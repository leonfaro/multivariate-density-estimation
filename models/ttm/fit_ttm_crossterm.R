## Cross-term Triangular Transport Map (TTM)
## Concrete fitter and evaluator; no proxies, R-only, deterministic.

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# Source required TTM components (relative paths)
local({
  p_b <- file.path("models", "ttm", "ttm_bases.R"); if (file.exists(p_b)) source(p_b)
  p_cb <- file.path("models", "ttm", "crossterm_basis.R"); if (file.exists(p_cb)) source(p_cb)
  p_c <- file.path("models", "ttm", "ttm_core.R");  if (file.exists(p_c)) source(p_c)
  p_o <- file.path("models", "ttm", "ttm_optimize.R"); if (file.exists(p_o)) source(p_o)
})

.std_stats <- function(X) {
  mu <- colMeans(X)
  sigma <- apply(X, 2, sd) + .Machine$double.eps
  list(mu = mu, sigma = sigma)
}

.standardize <- function(X, mu, sigma) {
  sweep(sweep(X, 2, mu, "-"), 2, sigma, "/")
}

#' Fit cross-term TTM via maps-from-samples with GL quadrature
#'
#' @param data matrix X or list with X_tr, X_te
#' @param seed integer RNG seed
#' @param deg_g degree for predecessor polynomial features g(x_<k>)
#' @param df_t B-spline df for h-basis in t=x_k
#' @param Q Gauss–Legendre nodes on [0,1]
#' @param lambda ridge weight (split internally into non/mon parts)
#' @param Hmax clamp for h before exponentiation
#' @param maxit optimizer max iterations (sets option cross.maxit)
#' @return list(S, NLL_train, NLL_test, time_train, time_pred)
fit_ttm_crossterm <- function(data, seed = 42L, deg_g = 2L, df_t = 6L, Q = 16L, lambda = 1e-3, Hmax = 8, maxit = NULL) {
  set.seed(as.integer(seed))
  if (!is.null(maxit)) options(cross.maxit = as.integer(maxit))
  # Accept split or matrix
  S_in <- if (is.list(data) && !is.null(data$X_tr) && !is.null(data$X_te)) data else {
    stopifnot(is.matrix(data))
    if (!exists("split_data")) source("02_split.R")
    list(X_tr = split_data(data, seed)$X_tr, X_te = split_data(data, seed)$X_te)
  }
  X_tr <- as.matrix(S_in$X_tr); X_te <- as.matrix(S_in$X_te)
  K <- ncol(X_tr); N <- nrow(X_tr)

  # Helper: oriented integral and clipping stats for given k
  .oriented_integral <- function(u, Uprev, beta, spec_h, nodes, weights, H, stats_env = NULL) {
    N <- length(u)
    I <- rep(0, N)
    max_abs_h <- -Inf
    max_exp_ht <- -Inf
    max_abs_h_t <- -Inf
    n_clip_pos <- 0L
    n_clip_neg <- 0L
    for (q in seq_along(nodes)) {
      tq <- u * nodes[q]
      Hq <- build_h(tq, Uprev, spec_h)
      hq <- as.numeric(Hq %*% beta)
      # Count clipping before exponentiation
      n_clip_pos <- n_clip_pos + sum(hq > H, na.rm = TRUE)
      n_clip_neg <- n_clip_neg + sum(hq < -H, na.rm = TRUE)
      max_abs_h <- max(max_abs_h, max(abs(hq), na.rm = TRUE))
      htil <- pmin(pmax(hq, -H), H)
      eh <- exp(htil)
      max_abs_h_t <- max(max_abs_h_t, max(abs(htil), na.rm = TRUE))
      max_exp_ht <- max(max_exp_ht, max(eh, na.rm = TRUE))
      I <- I + weights[q] * eh
    }
    if (!is.null(stats_env)) {
      stats_env$n_clip_pos <- stats_env$n_clip_pos + as.integer(n_clip_pos)
      stats_env$n_clip_neg <- stats_env$n_clip_neg + as.integer(n_clip_neg)
      stats_env$max_abs_h  <- max(stats_env$max_abs_h, max_abs_h)
      stats_env$max_exp_ht <- max(stats_env$max_exp_ht, max_exp_ht)
      stats_env$max_abs_h_t <- max(stats_env$max_abs_h_t, max_abs_h_t)
    }
    u * I
  }

  # Prepare symmetric t-basis spec from empirical t-samples across train
  gl <- gauss_legendre_nodes(as.integer(Q))
  nodes <- gl$nodes; weights <- gl$weights
  # Standardize using training to get U
  st0 <- .std_stats(as.matrix(if (is.list(data) && !is.null(data$X_tr)) data$X_tr else data))
  mu0 <- st0$mu; sigma0 <- st0$sigma
  X0 <- if (is.list(data) && !is.null(data$X_tr)) as.matrix(data$X_tr) else as.matrix(data)
  U0 <- .standardize(X0, mu0, sigma0)
  # Sample t = u * xi across dims and nodes (subsample rows for efficiency)
  set.seed(as.integer(seed))
  N0 <- nrow(U0); K0 <- ncol(U0)
  idx <- if (N0 > 200) sample.int(N0, 200) else seq_len(N0)
  U_s <- as.matrix(U0[idx, , drop = FALSE])
  t_list <- lapply(seq_len(K0), function(k) as.numeric(U_s[, k]) %o% nodes)
  t_sample <- as.numeric(do.call(c, lapply(t_list, as.vector)))
  # Build symmetric knots and boundary; choose effective df as option-capped
  df_opt <- .get_cross_df_t(df_t)
  df_eff <- as.integer(min(df_opt, max(4L, floor(nrow(U_s) / 10))))
  deg_bs <- 3L
  sk <- build_symmetric_knots(t_sample, df = df_eff, degree = deg_bs)
  # Effective g-degree (align with optimizer default cap)
  deg_g_eff <- max(as.integer(deg_g), .get_cross_deg_g())
  spec_h_global <- list(df = df_eff, degree = deg_bs, deg_g = as.integer(deg_g_eff),
                        knots_t = sk$interior, boundary_t = sk$boundary, tail_const = TRUE)

  # Train
  time_train <- system.time({
    st <- .std_stats(X_tr); mu <- st$mu; sigma <- st$sigma
    Xs <- .standardize(X_tr, mu, sigma)
    coeffs <- vector("list", K)
    spec_h <- spec_h_global; gl_nodes <- nodes; gl_weights <- weights
    # Clip event log accumulator
    clip_log <- list()
    for (k in seq_len(K)) {
      x_prev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
      xk <- Xs[, k]
      optk <- .opt_crossterm_k(x_prev, xk, deg_g = deg_g_eff, df_t = df_t, lambda = lambda, Q = Q, Hmax = Hmax, k = k,
                               spec_h = spec_h_global, nodes_override = nodes, weights_override = weights)
      coeffs[[k]] <- list(alpha = optk$alpha, beta = optk$beta)
      # spec_h already set globally
      # Post-fit clipping statistics on training set for this k
      stats_env <- new.env(parent = emptyenv())
      stats_env$n_clip_pos <- 0L
      stats_env$n_clip_neg <- 0L
      stats_env$max_abs_h  <- -Inf
      stats_env$max_exp_ht <- -Inf
      stats_env$max_abs_h_t <- -Inf
      invisible(.oriented_integral(xk, x_prev, coeffs[[k]]$beta, spec_h, gl_nodes, gl_weights, Hmax, stats_env))
      clip_log[[k]] <- data.frame(k = k,
                                  iter = NA_integer_,
                                  n_clipped_pos = as.integer(stats_env$n_clip_pos),
                                  n_clipped_neg = as.integer(stats_env$n_clip_neg),
                                  max_abs_h_raw = as.numeric(stats_env$max_abs_h),
                                  max_abs_h = as.numeric(stats_env$max_abs_h_t),
                                  max_exp_h_tilde = as.numeric(stats_env$max_exp_ht))
      msg <- sprintf("[CTM][k=%d] n_clipped_pos=%d n_clipped_neg=%d max|h|=%.6g max|h_raw|=%.6g max exp(h~)=%.6g",
                     k, clip_log[[k]]$n_clipped_pos, clip_log[[k]]$n_clipped_neg,
                     clip_log[[k]]$max_abs_h, clip_log[[k]]$max_abs_h_raw, clip_log[[k]]$max_exp_h_tilde)
      try(message(msg), silent = TRUE)
    }
    S <- list(
      algo = "crossterm",
      mu = mu, sigma = sigma,
      coeffs = coeffs,
      spec_h = spec_h,
      gl_nodes = gl_nodes,
      gl_weights = gl_weights,
      Hmax = Hmax,
      deg_g = as.integer(deg_g_eff)
    )
    class(S) <- "ttm_cross_term"
    assign(".last_fit_ttm_crossterm", S, envir = .GlobalEnv)
    # Persist clip log CSV
    dir.create("artifacts", showWarnings = FALSE)
    cldf <- do.call(rbind, clip_log)
    utils::write.csv(cldf, file = file.path("artifacts", "cross_clip_events.csv"), row.names = FALSE)
    # Basis summary logging
    # Summarize t-basis on sample t_sample
    Bt_s <- .bspline_basis_open_knots(t_sample, interior = spec_h$knots_t %||% numeric(0), degree = spec_h$degree %||% 3L, boundary = spec_h$boundary_t %||% range(t_sample))
    if (isTRUE(spec_h$tail_const %||% TRUE)) Bt_s <- cbind(Bt_s, tail_plateau_features(t_sample, spec_h$boundary_t))
    Gx_s <- if (ncol(Xs) > 1) build_g(matrix(0, nrow = length(t_sample), ncol = ncol(Xs) - 1), deg = spec_h$deg_g) else matrix(1, nrow = length(t_sample), ncol = 1)
    H_dim <- ncol(Bt_s) * ncol(Gx_s)
    basis_df <- data.frame(
      t_basis_dim = ncol(Bt_s), g_dim = ncol(Gx_s), h_dim = H_dim,
      degree = spec_h$degree, df = spec_h$df, boundary_left = spec_h$boundary_t[1], boundary_right = spec_h$boundary_t[2],
      stringsAsFactors = FALSE
    )
    # Add simple means/sds for t-basis columns
    cm <- colMeans(Bt_s); cs <- apply(Bt_s, 2, sd)
    # Write CSV with first few moments and knots
    kn <- spec_h$knots_t %||% numeric(0)
    basis_df$knots <- paste(sprintf("%.6g", kn), collapse = ",")
    basis_df$colmean_head <- paste(sprintf("%.6g", head(cm, min(5, length(cm)))), collapse = ",")
    basis_df$colsd_head <- paste(sprintf("%.6g", head(cs, min(5, length(cs)))), collapse = ",")
    utils::write.csv(basis_df, file = file.path("artifacts", "cross_basis_dims.csv"), row.names = FALSE)
  })[["elapsed"]]

  # Predict-once timing and numeric sanity
  time_pred <- system.time({ invisible(predict_ttm(S, X_te, type = "logdensity_by_dim")) })[["elapsed"]]
  LD_tr <- predict_ttm(S, X_tr, type = "logdensity_by_dim")
  LD_te <- predict_ttm_crossterm(S, X_te, type = "logdensity_by_dim")
  stopifnot(is.matrix(LD_tr), is.matrix(LD_te), all(is.finite(LD_tr)), all(is.finite(LD_te)))

  # Jacobian positivity check (train/test)
  J_tr <- predict_ttm(S, X_tr, type = "jac_diag")
  J_te <- predict_ttm(S, X_te, type = "jac_diag")
  if (any(!(J_tr > 0)) || any(!(J_te > 0))) {
    bad <- which(!(J_te > 0), arr.ind = TRUE)
    if (length(bad) == 0) bad <- which(!(J_tr > 0), arr.ind = TRUE)
    i <- bad[1, 1]; k <- bad[1, 2]
    # Recompute raw h for offending point on test
    mu <- S$mu; sigma <- S$sigma
    U_te <- .standardize(X_te, mu, sigma)
    u <- U_te[i, k]
    Uprev <- if (k > 1) U_te[i, 1:(k - 1), drop = FALSE] else matrix(0, 1, 0)
    Hstar <- build_h(u, Uprev, S$spec_h)
    h_val <- as.numeric(Hstar %*% S$coeffs[[k]]$beta)
    stop(sprintf("Nonpositive Jacobian detected at i=%d, k=%d, u=%.6g, h(u)=%.6g", i, k, u, h_val))
  }

  # Oriented sign checks on held-out set
  U_te <- .standardize(X_te, S$mu, S$sigma)
  N_te <- nrow(U_te)
  rows <- list()
  for (k in seq_len(ncol(U_te))) {
    Uprev <- if (k > 1) U_te[, 1:(k - 1), drop = FALSE] else matrix(0, N_te, 0)
    u <- U_te[, k]
    Iu <- .oriented_integral(u, Uprev, S$coeffs[[k]]$beta, S$spec_h, S$gl_nodes, S$gl_weights, S$Hmax)
    int_code <- sign(Iu)
    sgn_u <- sign(u)
    sign_ok <- (int_code == sgn_u) | ((u == 0) & (Iu == 0))
    rows[[k]] <- data.frame(i = seq_len(N_te), k = k, u_k = as.numeric(u), int_code = as.integer(int_code), sign_ok = as.logical(sign_ok))
  }
  sign_df <- do.call(rbind, rows)
  dir.create("artifacts", showWarnings = FALSE)
  utils::write.csv(sign_df, file = file.path("artifacts", "cross_sign_checks.csv"), row.names = FALSE)

  # Monotonicity-of-integral check: u grid, per-sample curves must be non-decreasing
  grid_u <- seq(-4, 4, length.out = 41)
  viol <- list()
  for (k in seq_len(ncol(U_te))) {
    Uprev <- if (k > 1) U_te[, 1:(k - 1), drop = FALSE] else matrix(0, N_te, 0)
    for (i in seq_len(N_te)) {
      up <- if (ncol(Uprev) > 0) matrix(Uprev[i, , drop = FALSE], nrow = length(grid_u), ncol = ncol(Uprev), byrow = TRUE) else matrix(0, length(grid_u), 0)
      Igrid <- .oriented_integral(grid_u, up, S$coeffs[[k]]$beta, S$spec_h, S$gl_nodes, S$gl_weights, S$Hmax)
      ok <- all(diff(Igrid) >= -1e-8)  # numerical tolerance
      if (!ok) viol[[length(viol) + 1L]] <- data.frame(i = i, k = k, violations = sum(diff(Igrid) < -1e-8), total_grid = length(grid_u))
    }
  }
  vio_df <- if (length(viol) > 0) do.call(rbind, viol) else data.frame(i = integer(0), k = integer(0), violations = integer(0), total_grid = integer(0))
  utils::write.csv(vio_df, file = file.path("artifacts", "cross_integral_checks.csv"), row.names = FALSE)
  if (nrow(vio_df) > 0) {
    msg <- sprintf("Monotonicity-of-integral violated in %d cases", nrow(vio_df))
    if (isTRUE(getOption("cross.strict_monotone", TRUE))) stop(msg) else {
      warning(msg)
      try(message(paste("[WARN]", msg)), silent = TRUE)
    }
  }

  # Tail slope report for |u|>4: slope b = exp(h) in tails
  for (k in seq_len(ncol(U_te))) {
    u <- U_te[, k]
    idx_tail <- which(abs(u) >= 4)
    if (length(idx_tail) > 0) {
      Uprev <- if (k > 1) U_te[idx_tail, 1:(k - 1), drop = FALSE] else matrix(0, length(idx_tail), 0)
      Hst <- build_h(u[idx_tail], Uprev, S$spec_h)
      h_raw <- as.numeric(Hst %*% S$coeffs[[k]]$beta)
      h_til <- pmin(pmax(h_raw, -S$Hmax), S$Hmax)
      b <- exp(h_til)
      try(message(sprintf("[TAIL][k=%d] |u|>=4 slope b range: [%.6g, %.6g]", k, min(b), max(b))), silent = TRUE)
    } else {
      try(message(sprintf("[TAIL][k=%d] no samples with |u|>=4 in test set", k)), silent = TRUE)
    }
  }

  list(S = S,
       NLL_train = mean(-predict_ttm(S, X_tr, type = "logdensity")),
       NLL_test  = mean(-predict_ttm_crossterm(S, X_te,  type = "logdensity")),
       time_train = time_train, time_pred = time_pred)
}

# Convenience evaluator; delegates to generic predict_ttm
predict_ttm_crossterm <- function(object, newdata, type = c("logdensity_by_dim", "logdensity")) {
  predict_ttm(object, newdata, match.arg(type))
}

# S3 predict wrapper for cross-term models
predict.ttm_cross_term <- function(object, newdata, type = c("logdensity_by_dim", "logdensity"), ...) {
  predict_ttm(object, newdata, match.arg(type))
}

# Compatibility wrapper mirroring separable/marginal helpers
trainCrossTermMap <- function(S, seed = 42, ...) {
  fit_ttm(S, algo = "crossterm", seed = seed, ...)
}

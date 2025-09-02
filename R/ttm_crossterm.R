# Cross-term triangular map: S_k(x) = g_k(x_<k>) + ∫_0^{x_k} exp(h_k(t, x_<k})) dt
# Minimal maps-from-samples fit using unified bases and core predictor.

if (!exists("build_f") || !exists("d_build_f") || !exists("build_g")) {
  if (exists("root_path")) {
    src1 <- file.path(root_path, "models", "ttm", "ttm_bases.R")
  } else {
    src1 <- file.path("models", "ttm", "ttm_bases.R")
  }
  if (file.exists(src1)) source(src1)
}
if (!exists("ttm_forward")) {
  if (exists("root_path")) {
    src2 <- file.path(root_path, "models", "ttm", "ttm_core.R")
  } else {
    src2 <- file.path("models", "ttm", "ttm_core.R")
  }
  if (file.exists(src2)) source(src2)
}

.std_stats <- function(X) {
  mu <- colMeans(X)
  sigma <- apply(X, 2, sd) + .Machine$double.eps
  list(mu = mu, sigma = sigma)
}

.standardize <- function(X, mu, sigma) sweep(sweep(X, 2, mu, "-"), 2, sigma, "/")

fit_ttm_crossterm <- function(data, deg_g = 2L, df_t = 6L, lambda = 1e-3, Q = 16L, Hmax = 20, maxit = 50L, seed = 42, ...) {
  set.seed(seed)
  suppressWarnings(try({ RNGkind("L'Ecuyer-CMRG") }, silent = TRUE))
  S_in <- if (is.list(data) && !is.null(data$X_tr) && !is.null(data$X_te)) data else {
    stopifnot(is.matrix(data))
    if (!exists("split_data")) source("02_split.R")
    list(X_tr = split_data(data, seed)$X_tr, X_te = split_data(data, seed)$X_te)
  }
  X_tr <- S_in$X_tr; X_te <- S_in$X_te
  N <- nrow(X_tr); K <- ncol(X_tr)

  time_train <- system.time({
    st <- .std_stats(X_tr); mu <- st$mu; sigma <- st$sigma
    Xs <- .standardize(X_tr, mu, sigma)
    Q_used <- .get_cross_quad_nodes(Q)
    gl <- gauss_legendre_nodes(Q_used)
    nodes <- gl$nodes; weights <- gl$weights
    # Choose t-spline df via option with N-aware cap; degree remains 3L
    df_opt <- .get_cross_df_t(df_t)
    df_eff <- as.integer(min(df_opt, max(4L, floor(N / 10))))
    df_t_alt <- as.integer(df_t)
    # h-basis degree in predecessors remains as provided (do not change here)
    spec_h <- list(df = df_eff, degree = 3L, deg_g = as.integer(deg_g))
    coeffs <- vector("list", K)
    init_coeffs <- vector("list", K)
    local_maxit <- .get_cross_maxit()
    # Optional warm-start: initialize from separable fit
    if (isTRUE(getOption("cross.warmstart_from_separable", FALSE))) {
      sep_deg <- getOption("cross.sep_degree_g", 2L)
      sep_lam <- getOption("cross.sep_lambda", 1e-3)
      sep_fit <- fit_ttm(list(X_tr = X_tr, X_te = X_te), algo = "separable", degree_g = as.integer(sep_deg), lambda = sep_lam, seed = seed)
      for (k in seq_len(K)) {
        Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
        xk <- Xs[, k]
        deg_g_used <- max(as.integer(deg_g), .get_cross_deg_g())
        Gx <- if (ncol(Xprev) > 0) build_g(Xprev, deg = deg_g_used) else matrix(1, N, 1L)
        m_g <- ncol(Gx)
        # t-only h-basis size per k
        m_h <- ncol(build_h(xk, Xprev, spec_h))
        # g-init via LS
        Pnon_sep <- if (k > 1) build_g(Xprev, deg = sep_fit$S$degree_g) else matrix(1, N, 1L)
        c_non_sep <- sep_fit$S$coeffs[[k]]$c_non
        g_sep_val <- if (ncol(Pnon_sep) > 0) as.numeric(Pnon_sep %*% c_non_sep) else rep(0, N)
        alpha0 <- tryCatch(solve(crossprod(Gx) + 1e-6 * diag(m_g), crossprod(Gx, g_sep_val)), error = function(e) rep(0, m_g))
        # h-init via LS on t-only
        c_mon_sep <- sep_fit$S$coeffs[[k]]$c_mon
        gl_local <- gauss_legendre_nodes(length(nodes))
        T <- as.numeric(outer(xk, gl_local$nodes, `*`))
        fprime <- as.numeric(d_build_f(T) %*% c_mon_sep)
        y <- log(pmax(1e-6, fprime))
        Bt <- .bspline_basis_uniform(T, df = spec_h$df, degree = spec_h$degree)
        beta_t <- tryCatch(solve(crossprod(Bt) + 1e-6 * diag(ncol(Bt)), crossprod(Bt, y)), error = function(e) rep(0, ncol(Bt)))
        beta0 <- rep(0, m_h)
        m_gx <- m_g
        for (j in seq_len(spec_h$df)) {
          idx <- (j - 1L) * m_gx + 1L
          if (idx <= length(beta0)) beta0[idx] <- beta_t[j]
        }
        init_coeffs[[k]] <- list(alpha = alpha0, beta = beta0)
      }
      # Evaluate initial cross NLL vs separable
      initS <- list(algo = "crossterm", mu = mu, sigma = sigma, coeffs = init_coeffs,
                    deg_g = as.integer(max(as.integer(deg_g), .get_cross_deg_g())), spec_h = spec_h,
                    gl_nodes = nodes, gl_weights = weights, Hmax = Hmax,
                    quad = list(Q = length(nodes), df_t = spec_h$df))
      class(initS) <- "ttm_cross_term"
      nll_sep <- mean(-predict_ttm(sep_fit$S, X_te, type = "logdensity"))
      nll_init <- mean(-predict_ttm(initS,    X_te, type = "logdensity"))
      if (isTRUE(getOption("cross.verbose", FALSE)) || interactive()) {
        message(sprintf("[CTM] warm-start: NLL_test sep=%.6f init-cross=%.6f", nll_sep, nll_init))
      }
    }

    train_one_k <- function(k) {
      Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
      xk <- Xs[, k]
      # Expand g-degree to at least cross.deg_g, but keep h-basis unchanged
      deg_g_used <- max(as.integer(deg_g), .get_cross_deg_g())
      Gx <- if (ncol(Xprev) > 0) build_g(Xprev, deg = deg_g_used) else matrix(1, N, 1L)
      m_g <- ncol(Gx)
      # Old feature count (no extra build): 1 + p*deg_g if p>0 else 1
      p_prev <- ncol(Xprev)
      m_g_alt <- if (p_prev > 0) 1L + p_prev * as.integer(deg_g) else 1L
      # h-basis uses the same poly degree in predecessors
      Hstar <- build_h(xk, Xprev, spec_h)  # N x (df * m_g_h)
      m_h <- ncol(Hstar)
      # Separate ridge weights for g/h
      ln <- .get_cross_lambda_non(Q_used, N, fallback = if (is.finite(lambda)) lambda else NA_real_)
      lm <- .get_cross_lambda_mon(ln, Q_used, N, fallback = if (is.finite(lambda)) lambda else NA_real_)
      # Objective using numeric gradient (L-BFGS-B on fn only)
      integ_calls <- 0L; integ_time <- 0
      fn_k <- function(theta) {
        a <- if (m_g > 0) theta[seq_len(m_g)] else numeric(0)
        b <- theta[(m_g + 1L):length(theta)]
        g <- if (m_g > 0) as.numeric(Gx %*% a) else rep(0, N)
        # Integral via GL on [0, xk]
        t0 <- proc.time()[3]
        I <- rep(0, N)
        for (q in seq_along(nodes)) {
          tq <- xk * nodes[q]
          Hq <- build_h(tq, Xprev, spec_h)
          v <- as.numeric(Hq %*% b)
          v <- pmax(pmin(v, Hmax), -Hmax)
          I <- I + weights[q] * exp(v)
        }
        I <- sign(xk) * abs(xk) * I
        integ_time <<- integ_time + (proc.time()[3] - t0)
        integ_calls <<- integ_calls + 1L
        S <- g + I
        h_star <- as.numeric(Hstar %*% b)
        h_star <- pmax(pmin(h_star, Hmax), -Hmax)
        0.5 * sum(S^2) - sum(h_star) + 0.5 * (ln * sum(a^2) + lm * sum(b^2))
      }
      grad_k <- NULL
      if (.use_cross_analytic_grad()) {
        grad_k <- function(theta) {
          a <- if (m_g > 0) theta[seq_len(m_g)] else numeric(0)
          b <- theta[(m_g + 1L):length(theta)]
          g <- if (m_g > 0) as.numeric(Gx %*% a) else rep(0, N)
          # Recompute I and S using same quadrature
          I <- rep(0, N)
          Hq_list <- vector("list", length(nodes))
          expv_list <- vector("list", length(nodes))
          for (q in seq_along(nodes)) {
            tq <- xk * nodes[q]
            Hq <- build_h(tq, Xprev, spec_h)
            v <- as.numeric(Hq %*% b)
            v <- pmax(pmin(v, Hmax), -Hmax)
            ev <- exp(v)
            I <- I + weights[q] * ev
            Hq_list[[q]] <- Hq
            expv_list[[q]] <- ev
          }
          I <- sign(xk) * abs(xk) * I
          S <- g + I
          # grad w.r.t a
          ga <- if (m_g > 0) as.numeric(t(Gx) %*% S) + ln * a else numeric(0)
          # grad w.r.t b: sum_i S_i * ∂I_i/∂b - sum_i Hstar_i + lm b
          gb <- rep(0, m_h)
          sgnx <- sign(xk) * abs(xk)
          for (q in seq_along(nodes)) {
            wfac <- weights[q] * expv_list[[q]] * sgnx * S
            gb <- gb + as.numeric(t(Hq_list[[q]]) %*% wfac)
          }
          gb <- gb - colSums(Hstar) + lm * b
          gtot <- c(ga, gb)
          if (any(!is.finite(gtot))) stop(sprintf("[CTM] NaN/Inf gradient at k=%d", k))
          # Optional gradient check on small subset
          if (isTRUE(getOption("cross.grad_check", FALSE))) {
            idx <- seq_len(min(2L, length(gtot)))
            eps <- 1e-6
            ng <- rep(0, length(idx))
            for (j in seq_along(idx)) {
              e <- rep(0, length(gtot)); e[idx[j]] <- eps
              ng[j] <- (fn_k(theta + e) - fn_k(theta - e)) / (2 * eps)
            }
            re <- max(abs(ng - gtot[idx]) / (abs(ng) + 1e-8))
            if (re > 1e-3) stop(sprintf("[CTM] grad-check failed at k=%d relErr=%.3e", k, re))
          }
          gtot
        }
      }
      theta0 <- if (!is.null(init_coeffs[[k]])) c(init_coeffs[[k]]$alpha, init_coeffs[[k]]$beta) else c(rep(0, m_g), rep(0, m_h))
      ctrl0 <- tryCatch(control, error = function(e) NULL)
      mx <- tryCatch(ctrl0$maxit, error = function(e) NULL) %||% .get_cross_maxit()
      control <- modifyList(ctrl0 %||% list(), list(maxit = mx))
      control$trace <- 0L
      t_opt0 <- proc.time()[3]
      opt <- if (is.null(grad_k)) optim(theta0, fn = fn_k, method = "L-BFGS-B", control = control)
             else optim(theta0, fn = fn_k, gr = grad_k, method = "L-BFGS-B", control = control)
      t_opt <- round((proc.time()[3] - t_opt0) * 1000)
      a_hat <- if (m_g > 0) opt$par[seq_len(m_g)] else numeric(0)
      b_hat <- opt$par[(m_g + 1L):length(opt$par)]
      if (isTRUE(getOption("cross.verbose", FALSE)) || interactive()) {
        # g share, ridge share, norms
        num <- sum(a_hat^2); den <- num + sum(b_hat^2)
        g_share <- if (den > 0) num / den else NA_real_
        ms <- round(integ_time * 1000)
        reg <- 0.5 * (ln * sum(a_hat^2) + lm * sum(b_hat^2))
        reg_share <- if (opt$value != 0) reg / opt$value else NA_real_
        message(sprintf("[CTM] k=%s end_loss=%.6f conv=%d Q=%d df_t_old=%d df_t_new=%d h_coefs=%d integ_ms=%d calls=%d opt_ms=%d f=%s g=%s g_feats_old=%d g_feats_new=%d g_share=%.3f lam_non=%.3g lam_mon=%.3g reg_share=%.3f ||g||2=%.3f ||h||2=%.3f",
                        k %||% "?", opt$value, opt$convergence, length(nodes), df_t_alt, spec_h$df, m_h, ms, integ_calls, t_opt,
                        paste0(opt$counts[["function"]]), paste0(opt$counts[["gradient"]]),
                        m_g_alt, m_g, g_share, ln, lm, reg_share, sum(a_hat^2), sum(b_hat^2)))
        # Quick probe: delta I with Q and Q+8 over 16 random points
        set.seed(123)
        idx <- if (N > 16) sample.int(N, 16L) else seq_len(N)
        b <- b_hat
        # I with current Q
        Iq <- rep(0, length(idx))
        for (q in seq_along(nodes)) {
          tq <- xk[idx] * nodes[q]
          Hq <- build_h(tq, Xprev[idx, , drop = FALSE], spec_h)
          v <- as.numeric(Hq %*% b)
          v <- pmax(pmin(v, Hmax), -Hmax)
          Iq <- Iq + weights[q] * exp(v)
        }
        Iq <- sign(xk[idx]) * abs(xk[idx]) * Iq
        # I with Q+8
        Qp <- length(nodes) + 8L
        glp <- gauss_legendre_nodes(Qp)
        Ip <- rep(0, length(idx))
        for (q in seq_along(glp$nodes)) {
          tq <- xk[idx] * glp$nodes[q]
          Hq <- build_h(tq, Xprev[idx, , drop = FALSE], spec_h)
          v <- as.numeric(Hq %*% b)
          v <- pmax(pmin(v, Hmax), -Hmax)
          Ip <- Ip + glp$weights[q] * exp(v)
        }
        Ip <- sign(xk[idx]) * abs(xk[idx]) * Ip
        d <- mean(abs(Iq - Ip))
        message(sprintf("[CTM] k=%s probe_dI(Q,%d)≈%.3e on 16 pts", k %||% "?", Qp, d))
      }
      list(alpha = a_hat, beta = b_hat,
           opt_info = list(convergence = opt$convergence, counts = opt$counts,
                           value = opt$value, maxit = control$maxit))
    }

    workers <- .get_cross_workers()
    if (!is.finite(workers) || is.na(workers) || workers < 1L) workers <- 1L
    k_idx <- seq_len(K)
    t_run0 <- proc.time()[3]
    results <- NULL
    try_mc <- try({ parallel::mclapply(k_idx, train_one_k, mc.cores = workers, mc.set.seed = TRUE) }, silent = TRUE)
    if (!inherits(try_mc, "try-error") && is.list(try_mc)) {
      results <- try_mc
    } else if (is.finite(workers) && workers > 1L) {
      # PSOCK fallback
      cl <- parallel::makeCluster(workers, type = "PSOCK")
      on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
      # Load dependencies on workers
      parallel::clusterEvalQ(cl, { source("models/ttm/ttm_bases.R"); source("models/ttm/ttm_core.R") })
      parallel::clusterExport(cl, varlist = c("Xs","N","K","deg_g","spec_h","Hmax","Q_used","nodes","weights","lambda","train_one_k",".get_cross_lambda_non",".get_cross_lambda_mon",".use_cross_analytic_grad"), envir = environment())
      try({ parallel::clusterSetRNGStream(cl, if (is.null(seed)) 42L else seed) }, silent = TRUE)
      results <- parallel::parLapply(cl, k_idx, train_one_k)
    } else {
      results <- lapply(k_idx, train_one_k)
    }
    t_run <- round((proc.time()[3] - t_run0) * 1000)
    if (isTRUE(getOption("cross.verbose", FALSE)) || interactive()) {
      message(sprintf("[CTM] trained %d components with %d cores in %d ms", K, workers, t_run))
    }
    for (k in k_idx) coeffs[[k]] <- results[[k]]
    S <- list(algo = "crossterm", mu = mu, sigma = sigma, coeffs = coeffs,
              # store g-degree used for prediction; h-degree remains in spec_h
              deg_g = as.integer(max(as.integer(deg_g), .get_cross_deg_g())), spec_h = spec_h,
              gl_nodes = nodes, gl_weights = weights, Hmax = Hmax,
              quad = list(Q = length(nodes), df_t = spec_h$df))
    class(S) <- "ttm_cross_term"
  })[["elapsed"]]
  if (isTRUE(getOption("cross.verbose", FALSE)) || interactive()) {
    message(sprintf("[CTM] total_train_sec=%.3f", time_train))
  }

  time_pred <- system.time({ invisible(predict_ttm(S, X_te, type = "logdensity_by_dim")) })[["elapsed"]]
  list(S = S,
       NLL_train = mean(-predict_ttm(S, X_tr, type = "logdensity")),
       NLL_test  = mean(-predict_ttm(S, X_te,  type = "logdensity")),
       time_train = time_train, time_pred = time_pred)
}

# S3 predict wrapper for cross-term TTM
predict.ttm_cross_term <- function(object, newdata, type = c("logdensity_by_dim", "logdensity"), ...) {
  predict_ttm(object, newdata, match.arg(type))
}

# Compatibility wrapper for existing scripts/tests.
# Accepts grid/tuning arguments used in scripts even if not all are used internally.
# - degree_g: polynomial degree for predecessors in both g and h
# - df_t:     B-spline df for t-direction
# - Q:        Gauss–Legendre nodes
# - lambda_non, lambda_mon: merged into a single lambda (average) for this minimal implementation
# - warmstart_from_separable, sep_degree_g, sep_lambda: accepted for API compatibility (no effect here)
trainCrossTermMap <- function(S,
                              degree_g = 2L,
                              seed = 42,
                              df_t = 6L,
                              Q = 16L,
                              lambda_non = NA_real_,
                              lambda_mon = NA_real_,
                              lambda = NA_real_,
                              Hmax = 20,
                              maxit = 50L,
                              warmstart_from_separable = FALSE,
                              sep_degree_g = 2L,
                              sep_lambda = 1e-3,
                              ...) {
  # Map inputs to separate ridge weights; fallback to single lambda if provided; else heuristics in trainer
  lam_non_out <- if (is.finite(lambda_non)) lambda_non else if (is.finite(lambda)) lambda else NA_real_
  lam_mon_out <- if (is.finite(lambda_mon)) lambda_mon else if (is.finite(lambda)) lambda else NA_real_
  # Optional warm-start placeholder (no-op in minimal implementation)
  # if (isTRUE(warmstart_from_separable)) {
  #   invisible(fit_ttm(S, algo = "separable", degree_g = sep_degree_g, lambda = sep_lambda, seed = seed))
  # }
  fit_ttm(S, algo = "crossterm", deg_g = degree_g, df_t = df_t, Q = Q,
          lambda = NA_real_, lam_non = lam_non_out, lam_mon = lam_mon_out,
          Hmax = Hmax, maxit = maxit, seed = seed)
}

# Optional CLI helpers for integration with optparse/argparse
cross_add_cli_option <- function(parser) {
  if (requireNamespace("optparse", quietly = TRUE) && inherits(parser, "OptionParser")) {
    parser <- optparse::add_option(parser, c("--cross_maxit"), type = "integer",
                                   default = 200L, help = "Max iterations for Cross-Term optim (default: %default)")
    parser <- optparse::add_option(parser, c("--cross_deg_g"), type = "integer",
                                   default = 3L, help = "Degree for g(x_<k) basis in Cross-Term (default: %default)")
    parser <- optparse::add_option(parser, c("--cross_quad_nodes"), type = "integer",
                                   default = 32L, help = "Gauss–Legendre nodes for Cross-Term integral (default: %default)")
    parser <- optparse::add_option(parser, c("--cross_df_t"), type = "integer",
                                   default = 8L, help = "Degrees of freedom for t-spline in Cross-Term (default: %default)")
    parser <- optparse::add_option(parser, c("--cross_lambda_non"), type = "double",
                                   default = NA_real_, help = "Ridge for g-coefficients (default heuristic)")
    parser <- optparse::add_option(parser, c("--cross_lambda_mon"), type = "double",
                                   default = NA_real_, help = "Ridge for h-coefficients (default 0.5*lambda_non)")
    parser <- optparse::add_option(parser, c("--cross_use_analytic_grad"), type = "logical",
                                   default = TRUE, help = "Use analytic gradients for Cross-Term objective (default: %default)")
    parser <- optparse::add_option(parser, c("--cross_workers"), type = "integer",
                                   default = NA_integer_, help = "Parallel workers for Cross-Term per-component training")
    return(parser)
  }
  if (requireNamespace("argparse", quietly = TRUE) && inherits(parser, "ArgumentParser")) {
    parser$add_argument("--cross_maxit", type = "integer", default = 200L,
                        help = "Max iterations for Cross-Term optim (default: 200)")
    parser$add_argument("--cross_deg_g", type = "integer", default = 3L,
                        help = "Degree for g(x_<k) basis in Cross-Term (default: 3)")
    parser$add_argument("--cross_quad_nodes", type = "integer", default = 32L,
                        help = "Gauss–Legendre nodes for Cross-Term integral (default: 32)")
    parser$add_argument("--cross_df_t", type = "integer", default = 8L,
                        help = "Degrees of freedom for t-spline in Cross-Term (default: 8)")
    parser$add_argument("--cross_lambda_non", type = "double", default = NA,
                        help = "Ridge for g-coefficients (default heuristic)")
    parser$add_argument("--cross_lambda_mon", type = "double", default = NA,
                        help = "Ridge for h-coefficients (default 0.5*lambda_non)")
    parser$add_argument("--cross_use_analytic_grad", type = "logical", default = TRUE,
                        help = "Use analytic gradients for Cross-Term objective (default: TRUE)")
    parser$add_argument("--cross_workers", type = "integer", default = NA,
                        help = "Parallel workers for Cross-Term per-component training")
    return(parser)
  }
  parser
}

cross_apply_cli_option <- function(opts) {
  if (!is.null(opts) && !is.null(opts$cross_maxit)) {
    options(cross.maxit = as.integer(opts$cross_maxit))
  }
  if (!is.null(opts) && !is.null(opts$cross_deg_g)) {
    options(cross.deg_g = as.integer(opts$cross_deg_g))
  }
  if (!is.null(opts) && !is.null(opts$cross_quad_nodes)) {
    options(cross.quad_nodes = as.integer(opts$cross_quad_nodes))
  }
  if (!is.null(opts) && !is.null(opts$cross_df_t)) {
    options(cross.df_t = as.integer(opts$cross_df_t))
  }
  if (!is.null(opts) && !is.null(opts$cross_lambda_non)) {
    options(cross.lambda_non = as.numeric(opts$cross_lambda_non))
  }
  if (!is.null(opts) && !is.null(opts$cross_lambda_mon)) {
    options(cross.lambda_mon = as.numeric(opts$cross_lambda_mon))
  }
  if (!is.null(opts) && !is.null(opts$cross_use_analytic_grad)) {
    options(cross.use_analytic_grad = isTRUE(as.logical(opts$cross_use_analytic_grad)))
  }
  if (!is.null(opts) && !is.null(opts$cross_workers) && is.finite(as.numeric(opts$cross_workers))) {
    options(cross.workers = as.integer(opts$cross_workers))
  }
  invisible(NULL)
}

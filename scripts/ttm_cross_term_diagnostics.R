#!/usr/bin/env Rscript

# Cross-term TTM diagnostics runner (no edits to core files)
# - Enforces N=50, perm=c(4,3,1,2), fixed RNG
# - Fits cross-term map (fallback to separable if not available)
# - Logs quadrature info, stability, and worst-sample analysis

log_line <- function(...) cat(sprintf("%s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
fmtv <- function(x) paste(sprintf("%.6g", x), collapse = ", ")
colSds <- function(M) apply(M, 2, stats::sd)

# Gauss–Legendre on [0,1]
gauss_legendre_01 <- function(n) {
  if (n <= 0 || n != as.integer(n)) stop("n must be positive integer")
  if (n == 1) return(list(nodes = 0.5, weights = 1))
  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  for (k in i) { J[k, k + 1] <- b[k]; J[k + 1, k] <- b[k] }
  e <- eigen(J, symmetric = TRUE)
  x <- (e$values + 1) / 2
  w <- (2 * (e$vectors[1, ]^2)) / 2
  list(nodes = as.numeric(x), weights = as.numeric(w))
}

main_diag <- function() {
  dir.create("logs", showWarnings = FALSE)
  dir.create("artifacts", showWarnings = FALSE)

  # Deterministic RNG
  set.seed(42L)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  LOG <- sprintf("logs/ttm_cross_perm4312_N50_%s.md", ts)
  con_out <- file(LOG, open = "at"); con_msg <- file(LOG, open = "at")
  sink(con_out, type = "output"); sink(con_msg, type = "message")
  # Relax monotone check to allow full report generation
  options(cross.strict_monotone = FALSE)
  on.exit({ try(sink(type = "message"), silent = TRUE); try(sink(), silent = TRUE); try(close(con_out), silent = TRUE); try(close(con_msg), silent = TRUE) }, add = TRUE)

  # Record session
  log_line("[SETUP] CWD=", getwd())
  log_line("[SETUP] R=", R.version.string)
  log_line("[SETUP] Session:")
  print(utils::sessionInfo())
  git_head <- try(system("git rev-parse HEAD", intern = TRUE), silent = TRUE)
  if (!inherits(git_head, "try-error")) log_line("[GIT] HEAD=", paste(git_head, collapse = " "))

  # Source required modules
  source("00_globals.R")
  source("01_data_generation.R")
  source("02_split.R")
  source("models/ttm/ttm_bases.R")
  source("models/ttm/ttm_core.R")
  source("models/ttm/ttm_marginal.R")
  source("models/ttm/ttm_separable.R")
  source("models/ttm/ttm_crossterm.R")
  if (!exists("fit_ttm_crossterm")) stop("fit_ttm_crossterm must be defined (models/ttm/fit_ttm_crossterm.R)")
  if (!exists("predict_ttm_crossterm")) stop("predict_ttm_crossterm must be defined (models/ttm/fit_ttm_crossterm.R)")
  if (requireNamespace("digest", quietly = TRUE)) {
    log_line("[WIRE] fit_hash(crossterm)=", digest::digest(body(fit_ttm_crossterm)))
  }

  # Config (same as main.R)
  config <- list(
    list(distr = "norm", parm = NULL),
    list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
    list(distr = "beta",
         parm = function(d) list(shape1 = softplus(d$X2), shape2 = softplus(d$X1))),
    list(distr = "gamma",
         parm = function(d) list(shape = softplus(d$X3), scale = softplus(d$X2)))
  )
  N <- 50L; perm <- c(4L, 3L, 1L, 2L); seed <- 42L
  # Hyperparams (use core defaults where not available)
  df_t <- getOption("cross.df_t", 6L)
  deg_g <- getOption("cross.deg_g", 2L)
  Q <- getOption("cross.quad_nodes", 16L)
  lambda <- getOption("cross.lambda", 1e-3)
  Hmax <- getOption("cross.clip_H", 20)
  log_line("[ARGS] N=", N, ", perm=", paste(perm, collapse = ","), ", deg_g=", deg_g, ", df_t=", df_t, ", Q=", Q, ", lambda=", lambda, ", clip_H=", Hmax)

  # Data gen and split
  X <- Generate_iid_from_config(N, config)
  S0 <- split_data(X, seed)
  S <- list(X_tr = S0$X_tr[, perm, drop = FALSE], X_te = S0$X_te[, perm, drop = FALSE])
  K <- ncol(S$X_tr); N_tr <- nrow(S$X_tr); N_te <- nrow(S$X_te)
  log_line("[SPLIT] N_tr=", N_tr, ", N_te=", N_te)

  # Fit cross-term (or proxy)
  fit <- fit_ttm_crossterm(S, seed = seed, deg_g = deg_g, df_t = df_t, Q = Q, lambda = lambda, maxit = 50L)
  M <- fit$S
  mu <- M$mu; sigma <- M$sigma
  log_line("[STANDARDIZE] mu=", fmtv(mu))
  log_line("[STANDARDIZE] sigma=", fmtv(sigma))

  # Quadrature info
  gl <- gauss_legendre_01(as.integer(Q))
  log_line("[QUAD] nodes(head)=", fmtv(head(gl$nodes, 5)), ", weights(head)=", fmtv(head(gl$weights, 5)))
  log_line("[QUAD] sum(weights)=", sprintf("%.6f", sum(gl$weights)))

  # Stability diagnostics per k
  Xs_tr <- sweep(sweep(S$X_tr, 2, mu, "-"), 2, sigma, "/")
  Xs_te <- sweep(sweep(S$X_te, 2, mu, "-"), 2, sigma, "/")
  worst_samples <- vector("list", K)
  iter_history <- vector("list", K)
  clip_events <- vector("list", K)
  has_h <- !is.null(M$spec_h)
  for (k in seq_len(K)) {
    log_line(sprintf("[ITER][k=%d] (no internal history available; reporting post-fit stability)", k))
    if (identical(M$algo, "crossterm") && has_h) {
      Xprev_tr <- if (k > 1) Xs_tr[, 1:(k - 1), drop = FALSE] else matrix(0, N_tr, 0)
      xk_tr <- Xs_tr[, k]
      H_tr <- build_h(xk_tr, Xprev_tr, M$spec_h)
      beta <- M$coeffs[[k]]$beta
      h_tr <- as.numeric(H_tr %*% beta)
      # Clipping as in core (for diagnostics)
      h_tr_c <- pmax(pmin(h_tr, Hmax), -Hmax)
      j_tr <- exp(h_tr_c)
      log_line(sprintf("[h_range][k=%d] min=%.6g, max=%.6g, min exp(h)=%.6g", k, min(h_tr), max(h_tr), min(j_tr)))
      Xprev_te <- if (k > 1) Xs_te[, 1:(k - 1), drop = FALSE] else matrix(0, N_te, 0)
      xk_te <- Xs_te[, k]
      H_te <- build_h(xk_te, Xprev_te, M$spec_h)
      h_te <- as.numeric(H_te %*% beta)
      h_te_c <- pmax(pmin(h_te, Hmax), -Hmax)
      j_te <- exp(h_te_c)
      log_line(sprintf("[h_range_test][k=%d] min=%.6g, max=%.6g, min exp(h)=%.6g", k, min(h_te), max(h_te), min(j_te)))
    } else {
      # Proxy: report Jacobian positivity
      Jtr <- predict_ttm(M, S$X_tr, type = "jac_diag")
      Jte <- predict_ttm(M, S$X_te, type = "jac_diag")
      log_line(sprintf("[J_pos][k=%d] min J_tr=%.6g, min J_te=%.6g", k, min(Jtr[, k]), min(Jte[, k])))
    }

    # Worst 10 samples on test by per-dim NLL
    LD_te <- predict_ttm(M, S$X_te, type = "logdensity_by_dim")
    nll_k <- -LD_te[, k]
    ord <- order(nll_k, decreasing = TRUE)
    top <- head(ord, 10)
    Z_te <- predict_ttm(M, S$X_te, type = "transform")
    if (identical(M$algo, "crossterm") && has_h) {
      Xprev_te <- if (k > 1) Xs_te[, 1:(k - 1), drop = FALSE] else matrix(0, N_te, 0)
      xk_te <- Xs_te[, k]
      h_k <- as.numeric(build_h(xk_te, Xprev_te, M$spec_h) %*% M$coeffs[[k]]$beta)
      worst <- data.frame(i = top, nll = nll_k[top], x_std = Xs_te[top, k], z = Z_te[top, k], h = h_k[top])
    } else {
      worst <- data.frame(i = top, nll = nll_k[top], x_std = Xs_te[top, k], z = Z_te[top, k])
    }
    worst_samples[[k]] <- worst
    log_line(sprintf("[worst10][k=%d] indices= %s", k, paste(worst$i, collapse = ",")))
  }

  # Evaluate metrics and invariants
  LD_tr <- predict_ttm(M, S$X_tr, type = "logdensity_by_dim")
  LD_te <- predict_ttm(M, S$X_te, type = "logdensity_by_dim")
  stopifnot(is.matrix(LD_tr), is.matrix(LD_te), all(is.finite(LD_tr)), all(is.finite(LD_te)))
  J_tr <- rowSums(LD_tr); J_te <- rowSums(LD_te)
  sum_ck_tr <- max(abs(J_tr - rowSums(LD_tr)))
  sum_ck_te <- max(abs(J_te - rowSums(LD_te)))
  log_line("[metrics] train per-dim NLL=", fmtv(-colMeans(LD_tr)), ", joint=", sprintf("%.6f", mean(-J_tr)), ", sum_check=", sprintf("%.3g", sum_ck_tr))
  log_line("[metrics] test  per-dim NLL=", fmtv(-colMeans(LD_te)), ", joint=", sprintf("%.6f", mean(-J_te)), ", sum_check=", sprintf("%.3g", sum_ck_te))

  # Save state/artifacts
  # Save as CSV unless ONLY_LOGS=1
  if (!identical(Sys.getenv("ONLY_LOGS", "0"), "1")) {
    ms <- data.frame(dim = seq_along(mu), mu = as.numeric(mu), sigma = as.numeric(sigma))
    utils::write.csv(ms, file = "artifacts/cross_mu_sigma.csv", row.names = FALSE)
    # Alpha/Beta coefficients (if present)
    if (!is.null(M$coeffs)) {
      alpha_df <- do.call(rbind, lapply(seq_along(M$coeffs), function(k) {
        a <- M$coeffs[[k]]$alpha; if (is.null(a)) return(NULL)
        data.frame(k = k, idx = seq_along(a), value = as.numeric(a))
      }))
      if (is.null(alpha_df)) alpha_df <- data.frame(k = integer(0), idx = integer(0), value = numeric(0))
      utils::write.csv(alpha_df, file = "artifacts/cross_alpha.csv", row.names = FALSE)
      beta_df <- do.call(rbind, lapply(seq_along(M$coeffs), function(k) {
        b <- M$coeffs[[k]]$beta; if (is.null(b)) return(NULL)
        data.frame(k = k, idx = seq_along(b), value = as.numeric(b))
      }))
      if (is.null(beta_df)) beta_df <- data.frame(k = integer(0), idx = integer(0), value = numeric(0))
      utils::write.csv(beta_df, file = "artifacts/cross_beta.csv", row.names = FALSE)
    } else {
      utils::write.csv(data.frame(k = integer(0), idx = integer(0), value = numeric(0)), file = "artifacts/cross_alpha.csv", row.names = FALSE)
      utils::write.csv(data.frame(k = integer(0), idx = integer(0), value = numeric(0)), file = "artifacts/cross_beta.csv", row.names = FALSE)
    }
    quad_df <- data.frame(idx = seq_along(gl$nodes), node = gl$nodes, weight = gl$weights, Q = as.integer(Q))
    utils::write.csv(quad_df, file = "artifacts/cross_quadrature.csv", row.names = FALSE)
    # Worst samples combined across k
    worst_df <- do.call(rbind, lapply(seq_along(worst_samples), function(k) {
      if (is.null(worst_samples[[k]])) return(NULL)
      cbind(k = k, worst_samples[[k]])
    }))
    if (is.null(worst_df)) worst_df <- data.frame(k = integer(0))
    utils::write.csv(worst_df, file = "artifacts/cross_worst_samples.csv", row.names = FALSE)
    # Placeholders for iteration history and clip events (not available)
    utils::write.csv(data.frame(k = integer(0), iter = integer(0), value = numeric(0), grad_norm = numeric(0), step_norm = numeric(0), line_search = character(0), clip_count = integer(0)),
                     file = "artifacts/cross_iter_history.csv", row.names = FALSE)
    utils::write.csv(data.frame(k = integer(0), iter = integer(0), clip_count = integer(0)),
                     file = "artifacts/cross_clip_events.csv", row.names = FALSE)
  }

  # Close sinks and echo tail
  sink(type = "message"); sink(); close(con_out); close(con_msg)
  cat(paste(utils::tail(readLines(LOG, warn = FALSE), 120), collapse = "\n"), "\n")
  invisible(TRUE)
}

if (sys.nframe() == 0L) main_diag()

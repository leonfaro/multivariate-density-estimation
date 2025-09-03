#!/usr/bin/env Rscript

# Separable TTM diagnostics runner (no core API changes)
# - Enforces N=50, perm=c(4,3,1,2), fixed RNG
# - Fits separable map only and logs analytic structures and metrics

log_line <- function(...) cat(sprintf("%s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
fmtv <- function(x) paste(sprintf("%.6g", x), collapse = ", ")
colSds <- function(M) apply(M, 2, stats::sd)
condnum <- function(A) {
  A <- as.matrix(A)
  if (any(!is.finite(A))) return(NA_real_)
  if (nrow(A) > 0 && ncol(A) > 0) {
    if (nrow(A) == ncol(A)) {
      # square: 2-norm condition number via svd ratio
      s <- try(svd(A, nu = 0, nv = 0)$d, silent = TRUE)
      if (inherits(s, "try-error")) return(NA_real_)
      if (length(s) == 0) return(NA_real_)
      return(max(s) / max(min(s), .Machine$double.eps))
    } else {
      # rectangular: ratio of largest to smallest nonzero singular value
      s <- try(svd(A, nu = 0, nv = 0)$d, silent = TRUE)
      if (inherits(s, "try-error")) return(NA_real_)
      s <- s[s > 0]
      if (length(s) == 0) return(NA_real_)
      return(max(s) / min(s))
    }
  }
  NA_real_
}

main_diag <- function() {
  dir.create("logs", showWarnings = FALSE)
  dir.create("artifacts", showWarnings = FALSE)

  # Deterministic RNG
  set.seed(42L)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  LOG <- sprintf("logs/ttm_separable_perm4312_N50_%s.log", ts)
  con_out <- file(LOG, open = "at")
  con_msg <- file(LOG, open = "at")
  sink(con_out, type = "output")
  sink(con_msg, type = "message")
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

  # Config (same as main.R)
  config <- list(
    list(distr = "norm", parm = NULL),
    list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
    list(distr = "beta",
         parm = function(d) list(shape1 = softplus(d$X2),
                                 shape2 = softplus(d$X1))),
    list(distr = "gamma",
         parm = function(d) list(shape = softplus(d$X3),
                                 scale = softplus(d$X2)))
  )
  N <- 50L
  perm <- c(4L, 3L, 1L, 2L)
  seed <- 42L
  deg_g <- 2L
  lambda <- 1e-3
  log_line("[ARGS] N=", N, ", perm=", paste(perm, collapse = ","), ", deg_g=", deg_g, ", lambda=", lambda)

  # Data gen and split
  X <- Generate_iid_from_config(N, config)
  S0 <- split_data(X, seed)
  S <- list(X_tr = S0$X_tr[, perm, drop = FALSE], X_te = S0$X_te[, perm, drop = FALSE])
  K <- ncol(S$X_tr); N_tr <- nrow(S$X_tr); N_te <- nrow(S$X_te)
  log_line("[SPLIT] N_tr=", N_tr, ", N_te=", N_te)

  # Fit separable
  fit <- fit_ttm(S, algo = "separable", degree_g = deg_g, lambda = lambda, seed = seed)
  M <- fit$S
  stopifnot(identical(M$algo, "separable"))
  mu <- M$mu; sigma <- M$sigma
  log_line("[STANDARDIZE] mu=", fmtv(mu))
  log_line("[STANDARDIZE] sigma=", fmtv(sigma))

  # Build designs per k on TRAIN (standardized)
  Xs_tr <- sweep(sweep(S$X_tr, 2, mu, "-"), 2, sigma, "/")
  Xs_te <- sweep(sweep(S$X_te, 2, mu, "-"), 2, sigma, "/")
  I_N <- diag(N_tr)

  P_non_dims <- integer(K)
  P_mon_dims <- integer(K)
  B_dims     <- integer(K)
  c_non_list <- vector("list", K)
  c_mon_list <- vector("list", K)
  diag_summ  <- vector("list", K)

  for (k in seq_len(K)) {
    x_prev_tr <- if (k > 1) Xs_tr[, 1:(k - 1), drop = FALSE] else matrix(0, N_tr, 0)
    xk_tr <- Xs_tr[, k]
    Phi_non <- if (ncol(x_prev_tr) > 0) build_g(x_prev_tr, deg = deg_g) else if (deg_g == 0L) matrix(1, N_tr, 1L) else matrix(0, N_tr, 0)
    Phi_mon <- build_f(xk_tr)
    B <- d_build_f(xk_tr)
    m_non <- ncol(Phi_non); m_mon <- ncol(Phi_mon)
    P_non_dims[k] <- m_non; P_mon_dims[k] <- m_mon; B_dims[k] <- ncol(B)
    # Orthogonalization structures
    if (m_non > 0) {
      M_k <- solve(crossprod(Phi_non) + lambda * diag(m_non), t(Phi_non))   # (m_non x N)
      A_k <- (I_N - Phi_non %*% M_k) %*% Phi_mon                            # (N x m_mon)
      D_k <- M_k %*% Phi_mon                                               # (m_non x m_mon)
    } else {
      M_k <- matrix(0, 0, N_tr)
      A_k <- Phi_mon
      D_k <- matrix(0, 0, m_mon)
    }

    # Optimization was done in fit; evaluate diagnostics at learned coeffs
    c_non <- M$coeffs[[k]]$c_non
    c_mon <- M$coeffs[[k]]$c_mon
    c_non_list[[k]] <- c_non
    c_mon_list[[k]] <- c_mon

    # f'(x) = (B %*% c_mon) > 0
    fp_tr <- as.numeric(B %*% c_mon)
    fp_te <- as.numeric(d_build_f(Xs_te[, k]) %*% c_mon)

    # g_k and f_k for summaries
    gk <- if (m_non > 0) as.numeric(Phi_non %*% c_non) else rep(0, N_tr)
    fk <- as.numeric(Phi_mon %*% c_mon)
    Zk_tr <- gk + fk

    diag_summ[[k]] <- list(
      P_non_dim = c(N_tr, m_non), P_mon_dim = c(N_tr, m_mon), B_dim = c(N_tr, ncol(B)),
      cond_Pnon = if (m_non > 0) condnum(crossprod(Phi_non) + lambda * diag(m_non)) else NA_real_,
      cond_Pmon = condnum(crossprod(Phi_mon)),
      cond_BtB  = condnum(crossprod(B)),
      cond_A    = condnum(A_k),
      cond_D    = if (m_non > 0) condnum(D_k) else NA_real_,
      min_fprime_tr = min(fp_tr), cnt_fprime_tr_le0 = sum(fp_tr <= 0), med_fprime_tr = stats::median(fp_tr),
      min_fprime_te = min(fp_te), cnt_fprime_te_le0 = sum(fp_te <= 0), med_fprime_te = stats::median(fp_te),
      med_Zk_tr = stats::median(Zk_tr), mean_Zk_tr = mean(Zk_tr)
    )

    # Log per-k block
    log_line(sprintf("[design_mats][k=%d] P_non=%dx%d, P_mon=%dx%d, B=%dx%d", k, N_tr, m_non, N_tr, m_mon, N_tr, ncol(B)))
    log_line(sprintf("[cond][k=%d] cond(P'P+lamI)=%.3g, cond(Pm'Pm)=%.3g, cond(B'B)=%.3g, cond(A)=%.3g, cond(D)=%.3g",
                     k, diag_summ[[k]]$cond_Pnon, diag_summ[[k]]$cond_Pmon, diag_summ[[k]]$cond_BtB, diag_summ[[k]]$cond_A, diag_summ[[k]]$cond_D))
    log_line(sprintf("[positivity][k=%d] min f' tr=%.6g (viol=%d), min f' te=%.6g (viol=%d)", k, diag_summ[[k]]$min_fprime_tr, diag_summ[[k]]$cnt_fprime_tr_le0, diag_summ[[k]]$min_fprime_te, diag_summ[[k]]$cnt_fprime_te_le0))
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

  # Z correlations and normality tests on TEST
  Z_te <- predict_ttm(M, S$X_te, type = "transform")
  C <- try(stats::cor(Z_te), silent = TRUE)
  if (!inherits(C, "try-error")) {
    log_line("[corr(Z)_test]")
    print(round(C, 3))
  }
  for (k in seq_len(K)) {
    p <- try(stats::shapiro.test(Z_te[, k])$p.value, silent = TRUE)
    if (!inherits(p, "try-error")) log_line(sprintf("[shapiro][k=%d] p=%.4g", k, p))
  }

  # Save state/artifacts
  # Save as CSV unless ONLY_LOGS=1
  if (!identical(Sys.getenv("ONLY_LOGS", "0"), "1")) {
    dims_df <- data.frame(k = seq_len(K), P_non = P_non_dims, P_mon = P_mon_dims, B = B_dims)
    utils::write.csv(dims_df, file = "artifacts/separable_fit_dims.csv", row.names = FALSE)
    cmon_df <- do.call(rbind, lapply(seq_len(K), function(k) {
      v <- as.numeric(c_mon_list[[k]])
      if (length(v) == 0) return(NULL)
      data.frame(k = k, idx = seq_along(v), value = v)
    }))
    if (is.null(cmon_df)) cmon_df <- data.frame(k = integer(0), idx = integer(0), value = numeric(0))
    utils::write.csv(cmon_df, file = "artifacts/separable_c_mon.csv", row.names = FALSE)
    cnon_df <- do.call(rbind, lapply(seq_len(K), function(k) {
      v <- as.numeric(c_non_list[[k]])
      if (length(v) == 0) return(NULL)
      data.frame(k = k, idx = seq_along(v), value = v)
    }))
    if (is.null(cnon_df)) cnon_df <- data.frame(k = integer(0), idx = integer(0), value = numeric(0))
    utils::write.csv(cnon_df, file = "artifacts/separable_c_non.csv", row.names = FALSE)
    # Metrics
    per_train <- data.frame(split = "train", k = seq_len(K), nll = as.numeric(-colMeans(LD_tr)))
    per_test  <- data.frame(split = "test",  k = seq_len(K), nll = as.numeric(-colMeans(LD_te)))
    utils::write.csv(rbind(per_train, per_test), file = "artifacts/separable_metrics_per_dim.csv", row.names = FALSE)
    joint_df <- data.frame(split = c("train","test"),
                           joint_mean_nll = c(mean(-J_tr), mean(-J_te)),
                           sum_check = c(sum_ck_tr, sum_ck_te))
    utils::write.csv(joint_df, file = "artifacts/separable_metrics_joint.csv", row.names = FALSE)
  }

  # Close sinks and echo tail
  sink(type = "message"); sink()
  close(con_out); close(con_msg)
  cat(paste(utils::tail(readLines(LOG, warn = FALSE), 100), collapse = "\n"), "\n")
  invisible(TRUE)
}

if (sys.nframe() == 0L) main_diag()

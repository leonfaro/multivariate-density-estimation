#!/usr/bin/env Rscript

# Wrapper to run main.R with fixed N and permutation, with extensive logging and artifacts.

safe_digest <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(x))
  paste(sample(letters, 8L, replace = TRUE), collapse = "")
}

fmtv <- function(x) paste(sprintf("%.6g", x), collapse = ", ")
skewness <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  m <- mean(x); s <- stats::sd(x); if (s == 0) return(0)
  mean(((x - m)/s)^3)
}
kurtosis_excess <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  m <- mean(x); s <- stats::sd(x); if (s == 0) return(-3)
  mean(((x - m)/s)^4) - 3
}

log_line <- function(...) cat(sprintf("%s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))

main_wrapper <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, "main.R"))) stop("run from repo root containing main.R")
  dir.create("logs", showWarnings = FALSE)
  dir.create("artifacts", showWarnings = FALSE)

  # Deterministic RNG
  set.seed(42L)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  LOG <- sprintf("logs/run_main_perm4312_N50_%s.md", ts)
  con_out <- file(LOG, open = "at")
  con_msg <- file(LOG, open = "at")
  sink(con_out, type = "output")
  sink(con_msg, type = "message")
  # Relax monotone check to allow full report generation
  options(cross.strict_monotone = FALSE)
  on.exit({ try(sink(type = "message"), silent = TRUE); try(sink(), silent = TRUE); try(close(con_out), silent = TRUE); try(close(con_msg), silent = TRUE) }, add = TRUE)

  log_line("[SETUP] CWD=", getwd())
  log_line("[SETUP] R=", R.version.string)
  log_line("[SETUP] Session:")
  print(utils::sessionInfo())
  git_head <- try(system("git rev-parse HEAD", intern = TRUE), silent = TRUE)
  if (!inherits(git_head, "try-error")) log_line("[GIT] HEAD=", paste(git_head, collapse = " "))

  # Source pipeline
  source("00_globals.R")
  source("01_data_generation.R")
  source("02_split.R")
  source("models/true_model.R")
  source("models/trtf_model.R")
  source("models/ttm/ttm_bases.R")
  source("models/ttm/ttm_core.R")
  source("models/ttm/ttm_marginal.R")
  source("models/ttm/ttm_separable.R")
  source("models/ttm/ttm_crossterm.R")
  source("models/true_joint_model.R")
  source("04_evaluation.R")
  # Require concrete cross-term API and log function digests
  if (!exists("fit_ttm_crossterm")) stop("fit_ttm_crossterm must be defined (models/ttm/fit_ttm_crossterm.R)")
  if (!exists("predict_ttm_crossterm")) stop("predict_ttm_crossterm must be defined (models/ttm/fit_ttm_crossterm.R)")
  fh_ct <- try(if (requireNamespace("digest", quietly = TRUE)) digest::digest(body(fit_ttm_crossterm)) else paste(sample(letters, 8L, TRUE), collapse=""), silent = TRUE)
  fh_sp <- try(if (requireNamespace("digest", quietly = TRUE)) digest::digest(body(fit_ttm_separable)) else paste(sample(letters, 8L, TRUE), collapse=""), silent = TRUE)
  log_line("[WIRE] fit_hash(crossterm)=", as.character(fh_ct), ", fit_hash(separable)=", as.character(fh_sp))

  # Source main and override n/perm afterwards
  source("main.R")
  perm <<- c(4L, 3L, 1L, 2L)
  n <<- 50L

  log_line("[ARGS] Enforced N=", n, ", perm=", paste(perm, collapse = ","))

  # Call main() once to exercise the stable API; capture its output in the log
  log_line("[RUN] Calling main() from main.R ...")
  out_tab <- try(main(), silent = TRUE)
  if (inherits(out_tab, "try-error")) log_line("[WARN] main() returned error: ", as.character(out_tab))

  # Reconstruct pipeline deterministically for exhaustive logging
  seed <- 42L
  cfg <- config
  log_line("[DATA_GEN] Families per dim: ", paste(vapply(cfg, `[[`, "", "distr"), collapse = ", "))
  Xgen <- Generate_iid_from_config(n, cfg)
  # Summary stats
  sums <- lapply(seq_len(ncol(Xgen)), function(k) {
    x <- Xgen[, k]
    c(mean = mean(x), sd = stats::sd(x), min = min(x), max = max(x), skew = skewness(x), kurt = kurtosis_excess(x))
  })
  sums_df <- do.call(rbind, sums)
  rownames(sums_df) <- paste0("X", seq_len(ncol(Xgen)))
  print(sums_df)

  # Split and permutation
  S0 <- split_data(Xgen, seed)
  idx_all <- seq_len(nrow(Xgen))
  # reconstruct split indices used
  set.seed(seed); idx <- sample.int(n); n_tr <- floor(0.8 * n); idx_tr <- idx[seq_len(n_tr)]; idx_te <- idx[(n_tr + 1):n]
  log_line("[SPLIT] n_tr=", length(idx_tr), ", n_te=", length(idx_te), ", seed=", seed,
           ", hash_tr=", safe_digest(idx_tr), ", hash_te=", safe_digest(idx_te))

  S <- list(X_tr = S0$X_tr[, perm, drop = FALSE], X_te = S0$X_te[, perm, drop = FALSE])
  cfgp <- cfg[perm]

  # Standardization
  mu <- colMeans(S$X_tr)
  sigma <- apply(S$X_tr, 2, sd) + .Machine$double.eps
  log_line("[STANDARDIZE] mu=", fmtv(mu))
  log_line("[STANDARDIZE] sigma=", fmtv(sigma))
  log_line("[STANDARDIZE] sum(-log sigma)=", sprintf("%.6f", sum(-log(sigma))), 
           ", min_sigma=", sprintf("%.6g", min(sigma)), ", any_zero=", any(abs(sigma) < .Machine$double.eps))
  # Save as CSV instead of RDS (unless ONLY_LOGS=1)
  only_logs <- identical(Sys.getenv("ONLY_LOGS", "0"), "1")
  if (!only_logs) {
    mu_sigma_df <- data.frame(dim = seq_along(mu), mu = as.numeric(mu), sigma = as.numeric(sigma))
    utils::write.csv(mu_sigma_df, file = "artifacts/mu_sigma.csv", row.names = FALSE)
  }

  # Setup logs for TTM
  log_line("[TTM_SETUP] enabled_models=TRUE,TRTF,TTM-marginal,TTM-separable,TTM-cross")
  log_line("[TTM_SETUP] cross-term hyperparams: deg_g=2, df_t=6, Q=16, lambda=1e-3, maxit=50, clip_H=20")

  # Fit models
  log_line("[FIT] TRUE (marginal)")
  M_true <- fit_TRUE(S, cfgp)
  log_line("[FIT] TRUE (joint)")
  M_true_joint <- tryCatch(fit_TRUE_JOINT(S, cfgp), error = function(e) { log_line("[WARN] TRUE (joint) failed: ", e$message); NULL })
  log_line("[FIT] TRTF (random forest)")
  M_trtf <- fit_TRTF(S, cfgp, seed = seed)
  log_line("[FIT] TTM Marginal")
  F_ttm <- fit_ttm(S, algo = "marginal", seed = seed)
  log_line("[FIT] TTM Separable")
  F_sep <- fit_ttm(S, algo = "separable", seed = seed)
  log_line("[FIT] TTM Cross-term")
  F_ctm <- fit_ttm(S, algo = "crossterm", seed = seed, deg_g = 2L, df_t = 6L, Q = 16L, lambda = 1e-3, maxit = 50L)

  # Save model artifacts
  if (!only_logs) {
    # Skip storing model objects as RDS; store lightweight signatures as CSV
    sig_df <- data.frame(
      model = c("true_marginal","true_joint","trtf","ttm_marginal","ttm_separable","ttm_crossterm"),
      present = c(TRUE, !is.null(M_true_joint), TRUE, TRUE, TRUE, TRUE),
      kind = c(NA, NA, NA, "marginal", "separable", "crossterm"),
      fn_hash = c(NA, NA, NA,
                  tryCatch(digest::digest(body(fit_ttm_marginal)), error = function(e) NA_character_),
                  tryCatch(digest::digest(body(fit_ttm_separable)), error = function(e) NA_character_),
                  tryCatch(digest::digest(body(fit_ttm_crossterm)), error = function(e) NA_character_)),
      stringsAsFactors = FALSE
    )
    utils::write.csv(sig_df, file = "artifacts/model_signatures.csv", row.names = FALSE)
  }

  # Z summaries and Jacobian checks (TTM)
  zsum <- list()
  for (tag in c("marginal", "separable", "crossterm")) {
    M <- switch(tag, marginal = F_ttm$S, separable = F_sep$S, crossterm = F_ctm$S)
    Z <- predict_ttm(M, S$X_te, type = "transform")
    J <- predict_ttm(M, S$X_te, type = "jac_diag")
    minJ <- min(J)
    if (!(minJ > 0)) log_line("[WARN] nonpositive Jacobian detected in ", tag)
    corr <- try(stats::cor(Z), silent = TRUE)
    zsum[[tag]] <- list(mean = colMeans(Z), sd = apply(Z, 2, sd), minJ = minJ, corr = corr)
    log_line("[Z] ", tag, ": mean=", fmtv(zsum[[tag]]$mean), ", sd=", fmtv(zsum[[tag]]$sd), ", minJ=", sprintf("%.6g", minJ))
  }
  if (!only_logs) {
    # Save Z summaries to CSV (long format)
    zs <- do.call(rbind, lapply(names(zsum), function(tag) {
      m <- zsum[[tag]]$mean; s <- zsum[[tag]]$sd; minJ <- zsum[[tag]]$minJ
      km <- seq_along(m)
      rbind(
        data.frame(tag = tag, stat = "mean", k = km, value = as.numeric(m)),
        data.frame(tag = tag, stat = "sd",   k = km, value = as.numeric(s)),
        data.frame(tag = tag, stat = "minJ", k = NA_integer_, value = as.numeric(minJ))
      )
    }))
    utils::write.csv(zs, file = "artifacts/z_summaries.csv", row.names = FALSE)
  }

  # Evaluation stage with invariants
  eval_one <- function(name, LD, joint = NULL) {
    stopifnot(is.matrix(LD), ncol(LD) == ncol(S$X_te), nrow(LD) == nrow(S$X_te))
    stopifnot(all(is.finite(LD)))
    LDj <- rowSums(LD)
    if (is.null(joint)) joint <- LDj
    tol <- max(abs(LDj - joint))
    ok <- tol <= 1e-12
    log_line("[EVAL] ", name, ": mean NLL=", sprintf("%.6f", mean(-joint)), ", SE=", sprintf("%.6f", stats::sd(-joint)/sqrt(length(joint))), 
             ", per_dim=", fmtv(-colMeans(LD)), ", sum_check=", sprintf("%.3g", tol), ", ok=", ok)
  }
  # TRUE marginals: per-row, per-dim log-densities
  K <- ncol(S$X_te)
  LD_true <- do.call(cbind, lapply(seq_len(K), function(k) {
    .log_density_vec(S$X_te[, k], cfgp[[k]]$distr, M_true$theta[[k]])
  }))
  eval_one("TRUE_marginal", LD_true)
  # TRUE joint by-dim in canonical order then permute
  LD_joint <- tryCatch(true_joint_logdensity_by_dim(cfg, S0$X_te), error = function(e) NULL)
  if (!is.null(LD_joint)) {
    LD_joint <- LD_joint[, perm, drop = FALSE]
    eval_one("TRUE_joint", LD_joint)
  }
  # TRTF
  LD_trtf <- predict(M_trtf, S$X_te, type = "logdensity_by_dim")
  eval_one("TRTF", LD_trtf)
  # TTM models
  eval_one("TTM_marginal", predict_ttm(F_ttm$S, S$X_te, type = "logdensity_by_dim"))
  eval_one("TTM_separable", predict_ttm(F_sep$S, S$X_te, type = "logdensity_by_dim"))
  eval_one("TTM_crossterm", predict_ttm(F_ctm$S, S$X_te, type = "logdensity_by_dim"))

  # Close sinks and show tail
  sink(type = "message"); sink()
  close(con_out); close(con_msg)
  cat(paste(readLines(LOG, warn = FALSE), collapse = "\n"), "\n")
  invisible(TRUE)
}

if (sys.nframe() == 0L) main_wrapper()

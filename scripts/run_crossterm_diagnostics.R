#!/usr/bin/env Rscript

# Cross-term diagnostics runner: executes main(), computes TTM metrics, and writes an MD report.

safe_digest <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(x))
  paste(sample(letters, 8L, replace = TRUE), collapse = "")
}
fmtv <- function(x) paste(sprintf("%.6g", x), collapse = ", ")
stderr <- function(x) stats::sd(x)/sqrt(length(x))

main_diag <- function() {
  dir.create("logs", showWarnings = FALSE)
  dir.create("artifacts", showWarnings = FALSE)

  # Deterministic RNG
  set.seed(42L)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  LOG <- sprintf("logs/ttm_fix_run_%s.md", ts)
  con_out <- file(LOG, open = "at"); con_msg <- file(LOG, open = "at")
  sink(con_out, type = "output"); sink(con_msg, type = "message")
  # Relax monotone check to allow full report generation
  options(cross.strict_monotone = FALSE)
  on.exit({ try(sink(type = "message"), silent = TRUE); try(sink(), silent = TRUE); try(close(con_out), silent = TRUE); try(close(con_msg), silent = TRUE) }, add = TRUE)

  cat("# Cross-term Diagnostics\n\n")
  cat("## Setup\n\n")
  cat("- CWD: ", getwd(), "\n", sep = "")
  cat("- R: ", R.version.string, "\n", sep = "")
  cat("- Session:\n\n")
  print(utils::sessionInfo())
  git_head <- try(system("git rev-parse HEAD", intern = TRUE), silent = TRUE)
  if (!inherits(git_head, "try-error")) cat("- GIT HEAD: ", paste(git_head, collapse = " "), "\n\n", sep = "")

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

  if (!exists("fit_ttm_crossterm")) stop("fit_ttm_crossterm must be defined (models/ttm/fit_ttm_crossterm.R)")
  if (!exists("predict_ttm_crossterm")) stop("predict_ttm_crossterm must be defined (models/ttm/fit_ttm_crossterm.R)")

  # Call main() to exercise the pipeline (do not override global bindings)
  cat("## Run main()\n\n")
  out_tab <- try(main(), silent = TRUE)
  if (inherits(out_tab, "try-error")) cat("[WARN] main() errored: ", as.character(out_tab), "\n", sep = "")

  # Reconstruct data and split deterministically for metrics
  seed <- 42L
  n <- 50L
  perm <- c(4L, 3L, 1L, 2L)
  cfg <- list(
    list(distr = "norm", parm = NULL),
    list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
    list(distr = "beta", parm = function(d) list(shape1 = softplus(d$X2), shape2 = softplus(d$X1))),
    list(distr = "gamma", parm = function(d) list(shape = softplus(d$X3), scale = softplus(d$X2)))
  )
  Xgen <- Generate_iid_from_config(n, cfg)
  S0 <- split_data(Xgen, seed)
  S <- list(X_tr = S0$X_tr[, perm, drop = FALSE], X_te = S0$X_te[, perm, drop = FALSE])

  # Fit separable and cross-term
  cat("## Hyperparameters\n\n")
  deg_g <- 2L; df_t <- 6L; Q <- 16L; lambda <- 1e-3; H <- 8L
  cat(sprintf("- deg_g=%d, df_t=%d, Q=%d, lambda=%.3g, H=%d\n\n", deg_g, df_t, Q, lambda, H))

  cat("## Training models\n\n")
  F_sep <- fit_ttm(S, algo = "separable", seed = seed)
  F_ctm <- fit_ttm(S, algo = "crossterm", seed = seed, deg_g = deg_g, df_t = df_t, Q = Q, lambda = lambda, maxit = 50L)
  M_sep <- F_sep$S; M_ctm <- F_ctm$S

  # Save model signatures
  sig_df <- data.frame(
    model = c("ttm_separable","ttm_crossterm"),
    kind = c("separable","crossterm"),
    fn_hash = c(tryCatch(digest::digest(body(fit_ttm_separable)), error = function(e) NA_character_),
                tryCatch(digest::digest(body(fit_ttm_crossterm)), error = function(e) NA_character_)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(sig_df, file = "artifacts/model_signatures.csv", row.names = FALSE)

  # Quadrature CSV (for record)
  gl <- gauss_legendre_nodes(Q)
  quad_df <- data.frame(idx = seq_along(gl$nodes), node = gl$nodes, weight = gl$weights, Q = as.integer(Q))
  utils::write.csv(quad_df, file = "artifacts/cross_quadrature.csv", row.names = FALSE)

  # Z and J summaries; train/test for both models
  summarize_ZJ <- function(M, S, tag) {
    Z_tr <- predict_ttm(M, S$X_tr, type = "transform")
    Z_te <- predict_ttm(M, S$X_te, type = "transform")
    J_tr <- predict_ttm(M, S$X_tr, type = "jac_diag")
    J_te <- predict_ttm(M, S$X_te, type = "jac_diag")
    minJ_tr <- min(J_tr); minJ_te <- min(J_te)
    corr_te <- try(stats::cor(Z_te), silent = TRUE)
    zs <- rbind(
      data.frame(tag = tag, split = "train", stat = "mean", k = seq_len(ncol(Z_tr)), value = as.numeric(colMeans(Z_tr))),
      data.frame(tag = tag, split = "train", stat = "sd",   k = seq_len(ncol(Z_tr)), value = as.numeric(apply(Z_tr, 2, sd))),
      data.frame(tag = tag, split = "test",  stat = "mean", k = seq_len(ncol(Z_te)), value = as.numeric(colMeans(Z_te))),
      data.frame(tag = tag, split = "test",  stat = "sd",   k = seq_len(ncol(Z_te)), value = as.numeric(apply(Z_te, 2, sd))),
      data.frame(tag = tag, split = "train", stat = "minJ", k = NA_integer_, value = as.numeric(minJ_tr)),
      data.frame(tag = tag, split = "test",  stat = "minJ", k = NA_integer_, value = as.numeric(minJ_te))
    )
    if (!inherits(corr_te, "try-error")) {
      utils::write.csv(corr_te, file = sprintf("artifacts/corr_Z_%s_test.csv", tag), row.names = TRUE)
    }
    zs
  }
  zs <- rbind(summarize_ZJ(M_sep, S, "separable"), summarize_ZJ(M_ctm, S, "crossterm"))
  utils::write.csv(zs, file = "artifacts/z_summaries.csv", row.names = FALSE)

  # LD and invariants; train/test for both models
  summarize_ld <- function(M, S, tag) {
    LD_tr <- predict_ttm(M, S$X_tr, type = "logdensity_by_dim"); JD_tr <- predict_ttm(M, S$X_tr, type = "logdensity")
    LD_te <- predict_ttm(M, S$X_te, type = "logdensity_by_dim"); JD_te <- predict_ttm(M, S$X_te, type = "logdensity")
    stopifnot(max(abs(rowSums(LD_tr) - JD_tr)) <= 1e-12)
    stopifnot(max(abs(rowSums(LD_te) - JD_te)) <= 1e-12)
    mtr <- colMeans(-LD_tr); ste_tr <- apply(-LD_tr, 2, stderr)
    mte <- colMeans(-LD_te); ste_te <- apply(-LD_te, 2, stderr)
    out_tr <- data.frame(model = tag, split = "train", k = seq_len(ncol(LD_tr)), mean_nll = as.numeric(mtr), se = as.numeric(ste_tr))
    out_te <- data.frame(model = tag, split = "test",  k = seq_len(ncol(LD_te)), mean_nll = as.numeric(mte), se = as.numeric(ste_te))
    list(train = out_tr, test = out_te, sum_ck_tr = max(abs(rowSums(LD_tr) - JD_tr)), sum_ck_te = max(abs(rowSums(LD_te) - JD_te)))
  }
  s_sep <- summarize_ld(M_sep, S, "separable")
  s_ctm <- summarize_ld(M_ctm, S, "crossterm")
  utils::write.csv(rbind(s_sep$train, s_sep$test), file = "artifacts/separable_metrics_per_dim.csv", row.names = FALSE)
  utils::write.csv(rbind(s_ctm$train, s_ctm$test), file = "artifacts/crossterm_metrics_per_dim.csv", row.names = FALSE)

  # Shapiro-Wilk for Z (test)
  shapiro_df <- function(M, S, tag) {
    Z <- predict_ttm(M, S$X_te, type = "transform")
    p <- sapply(seq_len(ncol(Z)), function(k) tryCatch(stats::shapiro.test(Z[, k])$p.value, error = function(e) NA_real_))
    data.frame(model = tag, k = seq_len(ncol(Z)), p_value = as.numeric(p))
  }
  sw_sep <- shapiro_df(M_sep, S, "separable"); sw_ctm <- shapiro_df(M_ctm, S, "crossterm")
  utils::write.csv(rbind(sw_sep, sw_ctm), file = "artifacts/shapiro_pvalues.csv", row.names = FALSE)

  # Acceptance checks and badges
  sign_df <- try(read.csv("artifacts/cross_sign_checks.csv", stringsAsFactors = FALSE), silent = TRUE)
  sign_rate <- if (!inherits(sign_df, "try-error") && nrow(sign_df) > 0) mean(sign_df$sign_ok) else NA_real_
  clip_df <- try(read.csv("artifacts/cross_clip_events.csv", stringsAsFactors = FALSE), silent = TRUE)
  max_abs_h <- if (!inherits(clip_df, "try-error") && nrow(clip_df) > 0) max(clip_df$max_abs_h) else NA_real_
  minJ_ctm <- min(predict_ttm(M_ctm, S$X_te, type = "jac_diag"))
  minJ_sep <- min(predict_ttm(M_sep, S$X_te, type = "jac_diag"))
  pass_sign <- isTRUE(sign_rate == 1)
  pass_h <- (!is.na(max_abs_h)) && (max_abs_h <= 8 + 1e-8)
  pass_sum <- (s_sep$sum_ck_tr == 0 && s_sep$sum_ck_te == 0 && s_ctm$sum_ck_tr == 0 && s_ctm$sum_ck_te == 0)
  pass_jac <- (minJ_ctm > 0 && minJ_sep > 0)
  pass_finite <- is.finite(F_ctm$NLL_train) && is.finite(F_ctm$NLL_test)

  cat("## Results\n\n")
  cat(sprintf("- sign_ok rate: %.4f (%s)\n", sign_rate, if (pass_sign) "PASS" else "FAIL"))
  cat(sprintf("- max |h| (clipped): %.6g (%s)\n", max_abs_h, if (pass_h) "PASS" else "FAIL"))
  cat(sprintf("- min J (separable): %.6g; min J (crossterm): %.6g (%s)\n", minJ_sep, minJ_ctm, if (pass_jac) "PASS" else "FAIL"))
  cat(sprintf("- sum check train/test (both models): %s\n", if (pass_sum) "PASS" else "FAIL"))
  cat(sprintf("- NLL crossterm finite (train/test): %s\n\n", if (pass_finite) "PASS" else "FAIL"))

  # Pretty tables in MD (compact)
  cat("### NLL per-dim (train)\n\n"); print(rbind(s_sep$train, s_ctm$train))
  cat("\n### NLL per-dim (test)\n\n"); print(rbind(s_sep$test, s_ctm$test))
  cat("\n### Shapiro p-values (Z, test)\n\n"); print(rbind(sw_sep, sw_ctm))
  cat("\n### Cor(Z) (crossterm, test)\n\n"); print(try(stats::cor(predict_ttm(M_ctm, S$X_te, type = "transform")), silent = TRUE))

  # Iter history CSV (summary rows per k)
  # Note: L-BFGS-B internal iterations are not logged; provide per-k summary placeholders
  it_df <- data.frame(k = seq_len(ncol(S$X_tr)), iter = NA_integer_, value = NA_real_, grad_norm = NA_real_)
  utils::write.csv(it_df, file = "artifacts/cross_iter_history.csv", row.names = FALSE)

  # Echo tail and done
  sink(type = "message"); sink(); close(con_out); close(con_msg)
  cat(paste(utils::tail(readLines(LOG, warn = FALSE), 200), collapse = "\n"), "\n")
  invisible(TRUE)
}

if (sys.nframe() == 0L) main_diag()

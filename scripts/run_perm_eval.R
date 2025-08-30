#!/usr/bin/env Rscript
# Run pipeline with a custom permutation and omit True (Joint) for non-default order

perm_arg <- commandArgs(trailingOnly = TRUE)
perm_vec <- if (length(perm_arg) >= 1L) as.integer(strsplit(perm_arg[[1]], ",", fixed = TRUE)[[1]]) else c(4L,1L,2L,3L)

source("main.R")

# Overwrite global perm used by main() logic
perm <- perm_vec

# Re-create the essential parts of main(), skipping True (Joint) which relies on the default triangular order
seed <- as.integer(Sys.getenv("SEED", 42))
set.seed(seed)
prep <- prepare_data(n, config, seed = seed)
S0 <- prep$S
S <- list(
  X_tr  = S0$X_tr[, perm, drop = FALSE],
  X_te  = S0$X_te[, perm, drop = FALSE]
)
cfg <- config[perm]

t_true_tr  <- system.time(mod_true      <- fit_TRUE(S, cfg))[['elapsed']]
t_joint_tr <- NA_real_
mod_true_joint <- NULL
t_trtf_tr  <- system.time(mod_trtf      <- fit_TRTF(S, cfg, seed = seed))[['elapsed']]
mod_ttm     <- trainMarginalMap(S, seed = seed);   t_ttm_tr <- mod_ttm$time_train
mod_ttm_sep <- trainSeparableMap(S, seed = seed);  t_sep_tr <- mod_ttm_sep$time_train
use_rev_nf <- Sys.getenv('MDE_CTM_USE_REV_NF','0') %in% c('1','true','TRUE')
if (use_rev_nf) {
  mod_ttm_cross <- trainCrossTermMapReverseNF(S, degree_g = 2, degree_t = 2, seed = seed)
  t_ct_tr <- NA_real_
} else {
  mod_ttm_cross <- trainCrossTermMap(S, degree_g = 3, seed = seed, warmstart_from_separable = TRUE); t_ct_tr <- mod_ttm_cross$time_train
}

if (!is.null(mod_ttm_cross$S$meta$ridge)) {
  rr <- mod_ttm_cross$S$meta$ridge
  cat(sprintf("[RIDGE] lambda_non=%.3g, lambda_mon=%.3g\n", rr[["lambda_non"]], rr[["lambda_mon"]]))
} else if (!is.null(mod_ttm_cross$S$meta$lambda_non) && !is.null(mod_ttm_cross$S$meta$lambda_mon)) {
  cat(sprintf("[RIDGE] lambda_non=%.3g, lambda_mon=%.3g\n", mod_ttm_cross$S$meta$lambda_non, mod_ttm_cross$S$meta$lambda_mon))
}

t_true_te  <- system.time(logL_TRUE(mod_true, S$X_te))[['elapsed']]
t_joint_te <- NA_real_
t_trtf_te  <- system.time(predict(mod_trtf, S$X_te, type = "logdensity_by_dim"))[['elapsed']]
t_ttm_te   <- mod_ttm$time_pred
t_sep_te   <- mod_ttm_sep$time_pred
t_ct_te    <- mod_ttm_cross$time_pred

mods <- list(
  true = mod_true,
  trtf = mod_trtf,
  ttm  = mod_ttm,
  ttm_sep = mod_ttm_sep,
  ttm_cross = mod_ttm_cross
)

# Local table builder without True (Joint)
stderr <- function(x) stats::sd(x)/sqrt(length(x))
calc_tab_no_true_joint <- function(models, config, X_te) {
  K <- length(config)
  ll_true <- matrix(NA_real_, nrow = nrow(X_te), ncol = K)
  for (k in seq_len(K)) {
    ll_vec <- .log_density_vec(X_te[, k], config[[k]]$distr, models$true$theta[[k]])
    ll_true[, k] <- -ll_vec
  }
  ll_trtf <- -predict(models$trtf, X_te, type = "logdensity_by_dim")
  mean_true <- colMeans(ll_true); se_true <- apply(ll_true, 2, stderr)
  mean_trtf <- colMeans(ll_trtf); se_trtf <- apply(ll_trtf, 2, stderr)

  if (!is.null(models$ttm)) {
    ll_ttm <- -predict(models$ttm$S, X_te, type = "logdensity_by_dim")
    mean_ttm <- colMeans(ll_ttm); se_ttm <- apply(ll_ttm, 2, stderr)
    total_nll_ttm <- rowSums(ll_ttm); se_sum_ttm <- stats::sd(total_nll_ttm)/sqrt(nrow(X_te))
  } else { mean_ttm <- rep(NA_real_, K); se_ttm <- rep(NA_real_, K); se_sum_ttm <- NA_real_ }
  if (!is.null(models$ttm_sep)) {
    ll_ttm_sep <- -predict(models$ttm_sep$S, X_te, type = "logdensity_by_dim")
    mean_sep <- colMeans(ll_ttm_sep); se_sep <- apply(ll_ttm_sep, 2, stderr)
    total_nll_sep <- rowSums(ll_ttm_sep); se_sum_sep <- stats::sd(total_nll_sep)/sqrt(nrow(X_te))
  } else { mean_sep <- rep(NA_real_, K); se_sep <- rep(NA_real_, K); se_sum_sep <- NA_real_ }
  if (!is.null(models$ttm_cross)) {
    ll_ttm_cross <- -predict(models$ttm_cross$S, X_te, type = "logdensity_by_dim")
    mean_cross <- colMeans(ll_ttm_cross); se_cross <- apply(ll_ttm_cross, 2, stderr)
    total_nll_cross <- rowSums(ll_ttm_cross); se_sum_cross <- stats::sd(total_nll_cross)/sqrt(nrow(X_te))
  } else { mean_cross <- rep(NA_real_, K); se_cross <- rep(NA_real_, K); se_sum_cross <- NA_real_ }

  total_nll_true <- rowSums(ll_true); se_sum_true <- stats::sd(total_nll_true)/sqrt(nrow(X_te))
  total_nll_trtf <- rowSums(ll_trtf); se_sum_trtf <- stats::sd(total_nll_trtf)/sqrt(nrow(X_te))
  fmt <- function(m, se) sprintf("%.2f ± %.2f", round(m, 2), round(2 * se, 2))
  tab <- data.frame(
    dim = as.character(seq_len(K)),
    distribution = sapply(config, `[[`, "distr"),
    `True (marginal)` = fmt(mean_true, se_true),
    `Random Forest` = fmt(mean_trtf, se_trtf),
    `Marginal Map` = fmt(mean_ttm, se_ttm),
    `Separable Map` = fmt(mean_sep, se_sep),
    `Cross-term Map` = fmt(mean_cross, se_cross),
    stringsAsFactors = FALSE
  )
  sum_row <- data.frame(
    dim = "k",
    distribution = "SUM",
    `True (marginal)` = fmt(sum(mean_true), se_sum_true),
    `Random Forest` = fmt(sum(mean_trtf), se_sum_trtf),
    `Marginal Map` = fmt(sum(mean_ttm),  se_sum_ttm),
    `Separable Map` = fmt(sum(mean_sep),  se_sum_sep),
    `Cross-term Map` = fmt(sum(mean_cross), se_sum_cross),
    stringsAsFactors = FALSE
  )
  rbind(tab, sum_row)
}

tab <- calc_tab_no_true_joint(mods, cfg, S$X_te)
cat(sprintf("n=%d\n", n))
print(tab)
cat(sprintf("Permutation order %s (train/test only)\n", paste(perm, collapse = ",")))
time_tab <- data.frame(
  model = c("True (marginal)", "Random Forest",
            "Marginal Map", "Separable Map", "Cross-term Map"),
  train_sec = c(t_true_tr, t_trtf_tr, t_ttm_tr, t_sep_tr, t_ct_tr),
  test_sec = c(t_true_te, t_trtf_te, t_ttm_te, t_sep_te, t_ct_te),
  stringsAsFactors = FALSE
)
time_tab$total_sec <- with(time_tab, train_sec + test_sec)
print(time_tab)

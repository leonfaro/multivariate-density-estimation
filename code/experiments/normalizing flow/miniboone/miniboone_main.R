#!/usr/bin/env Rscript

# MiniBooNE runner without artifacts/logs. Expects pre-split CSVs under
# `experiments/normalizing flow/miniboone/` by default.
Sys.setenv(NO_ARTIFACTS = "1")
options(mde.no_artifacts = TRUE)

if (!exists("initialize_repo")) {
  source(file.path(getwd(), "R", "loader.R"))
}
root_path <- initialize_repo()

stderr <- function(x) stats::sd(x) / sqrt(length(x))
fmt <- function(m, se) sprintf("%.2f ± %.2f", round(m, 2), round(2 * se, 2))

read_num_csv <- function(path) {
  if (!file.exists(path)) stop(sprintf("CSV not found: %s", path))
  df <- utils::read.csv(path, header = TRUE, check.names = TRUE)
  df <- df[vapply(df, is.numeric, logical(1))]
  as.matrix(df)
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  m <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(m)) sub(paste0("^--", key, "="), "", m[1]) else default
}

default_dir <- file.path(root_path, "experiments", "normalizing flow", "miniboone")
train_csv <- get_arg("train_csv", Sys.getenv("MINIBOONE_TRAIN", file.path(default_dir, "miniboone_train.csv")))
val_csv   <- get_arg("val_csv",   Sys.getenv("MINIBOONE_VAL",   file.path(default_dir, "miniboone_val.csv")))
test_csv  <- get_arg("test_csv",  Sys.getenv("MINIBOONE_TEST",  file.path(default_dir, "miniboone_test.csv")))
seed <- as.integer(get_arg("seed", Sys.getenv("SEED", "42")))

set.seed(seed)
Xtr_all <- read_num_csv(train_csv)
Xva_all <- read_num_csv(val_csv)
Xte_all <- read_num_csv(test_csv)
common <- Reduce(intersect, list(colnames(Xtr_all), colnames(Xva_all), colnames(Xte_all)))
if (!length(common)) stop("No common numeric columns across MiniBooNE splits")
X_tr <- Xtr_all[, common, drop = FALSE]
X_va <- Xva_all[, common, drop = FALSE]
X_te <- Xte_all[, common, drop = FALSE]

S <- list(X_tr = X_tr, X_te = X_te)
cfg <- replicate(ncol(X_tr), list(distr = "norm"), simplify = FALSE)

# Fit models (robust to TRTF absence)
mods <- list()
mods$ttm <- fit_ttm(S, algo = "marginal", seed = seed)$S
mods$ttm_sep <- fit_ttm(S, algo = "separable", seed = seed)$S
mods$trtf <- tryCatch(fit_TRTF(S, cfg, seed = seed), error = function(e) NULL)

N <- nrow(X_te); K <- ncol(X_te)
tab <- data.frame(dim = as.character(seq_len(K)), distribution = rep("data", K), stringsAsFactors = FALSE)

if (!is.null(mods$trtf)) {
  ll_trtf <- tryCatch(-predict(mods$trtf, X_te, type = "logdensity_by_dim"), error = function(e) matrix(NA_real_, N, K))
  tab[["Random Forest"]] <- fmt(colMeans(ll_trtf), apply(ll_trtf, 2, stderr))
}
ll_ttm <- -predict_ttm(mods$ttm, X_te, type = "logdensity_by_dim")
tab[["Marginal Map"]] <- fmt(colMeans(ll_ttm), apply(ll_ttm, 2, stderr))
ll_sep <- -predict_ttm(mods$ttm_sep, X_te, type = "logdensity_by_dim")
tab[["Separable Map"]] <- fmt(colMeans(ll_sep), apply(ll_sep, 2, stderr))

sum_row <- data.frame(dim = "k", distribution = "SUM", stringsAsFactors = FALSE)
if (!is.null(mods$trtf)) sum_row[["Random Forest"]] <- fmt(sum(colMeans(ll_trtf)), stats::sd(rowSums(ll_trtf)) / sqrt(N))
sum_row[["Marginal Map"]]  <- fmt(sum(colMeans(ll_ttm)), stats::sd(rowSums(ll_ttm)) / sqrt(N))
sum_row[["Separable Map"]] <- fmt(sum(colMeans(ll_sep)), stats::sd(rowSums(ll_sep)) / sqrt(N))
tab <- rbind(tab, sum_row)

cat(sprintf("MiniBooNE dims=%d | Train/Val/Test CSVs loaded\n", K))
print(tab)

#!/usr/bin/env Rscript

# Compute PIT-based KS distances for the 4D generator across seeds

# Load repo and experiment config
ensure_loader <- function() {
  cur <- getwd()
  cands <- c(file.path(cur, "R", "loader.R"), file.path(cur, "code", "R", "loader.R"),
             file.path(dirname(cur), "R", "loader.R"), file.path(dirname(cur), "code", "R", "loader.R"))
  cand <- cands[file.exists(cands)][1]
  if (!length(cand)) stop("Could not locate loader.R")
  source(cand)
  initialize_repo()
}
root_path <- ensure_loader()
src <- function(...) source(file.path(root_path, ...))
src("experiments", "4d", "main.R")
src("R", "evaluation.R")

# Seeds and N
seeds <- Sys.getenv("SEEDS", "7,11,19")
seeds <- as.integer(strsplit(seeds, ",", fixed = TRUE)[[1]])
seeds <- seeds[is.finite(seeds)]
if (!length(seeds)) seeds <- c(7L, 11L, 19L)
N <- as.integer(Sys.getenv("N_OVERRIDE", "250"))
if (!is.finite(N) || N <= 0) N <- 250L

results_dir <- file.path(root_path, "experiments", "4d", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

rows <- list()
for (sd in seeds) {
  set.seed(sd)
  prep <- prepare_data(N, config, seed = sd)
  S <- prep$S
  # Fit models
  M_true <- fit_TRUE(S, config)
  M_trtf <- tryCatch(fit_TRTF(S, config, seed = sd), error = function(e) NULL)
  fit_m <- fit_ttm(S, algo = "marginal", seed = sd)
  fit_s <- fit_ttm(S, algo = "separable", seed = sd)
  mods <- list(true = M_true, trtf = M_trtf, ttm = fit_m$S, ttm_sep = fit_s$S)
  ksdf <- calc_ks_summary(mods, S, config = config)
  ksdf$seed <- sd
  rows[[length(rows) + 1]] <- ksdf
}

per_seed <- do.call(rbind, rows)
if ("ks_per_dim" %in% names(per_seed)) per_seed$ks_per_dim <- NULL

agg <- do.call(rbind, lapply(split(per_seed, per_seed$model), function(df) {
  m <- mean(df$ks_median, na.rm = TRUE)
  se <- stats::sd(df$ks_median, na.rm = TRUE) / sqrt(sum(is.finite(df$ks_median)))
  data.frame(model = df$model[1], ks_median_mean = m, ks_median_se = se, stringsAsFactors = FALSE)
}))
row.names(agg) <- NULL

utils::write.csv(per_seed, file = file.path(results_dir, "ks_4d_per_seed.csv"), row.names = FALSE)
utils::write.csv(agg, file = file.path(results_dir, "ks_4d_summary.csv"), row.names = FALSE)

print(agg)

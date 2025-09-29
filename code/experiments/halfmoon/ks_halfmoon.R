#!/usr/bin/env Rscript

# Compute PIT-based KS distances for half-moon across seeds

loader_path <- if (exists("locate_repo_loader")) locate_repo_loader() else {
  # best-effort search
  cur <- getwd()
  cands <- c(file.path(cur, "R", "loader.R"), file.path(cur, "code", "R", "loader.R"),
             file.path(dirname(cur), "R", "loader.R"), file.path(dirname(cur), "code", "R", "loader.R"))
  cand <- cands[file.exists(cands)][1]
  if (!length(cand)) stop("Could not locate loader.R from working directory")
  cand
}
source(loader_path)
root_path <- initialize_repo()

src <- function(...) source(file.path(root_path, ...))
src("experiments", "halfmoon", "halfmoon_data.R")
src("experiments", "halfmoon", "plot.R")

seeds <- Sys.getenv("SEEDS", "7,11,19")
seeds <- as.integer(strsplit(seeds, ",", fixed = TRUE)[[1]])
seeds <- seeds[is.finite(seeds)]
if (!length(seeds)) seeds <- c(7L, 11L, 19L)

N <- as.integer(Sys.getenv("N_OVERRIDE", "250"))
if (!is.finite(N) || N <= 0) N <- 250L
NOISE <- as.numeric(Sys.getenv("NOISE_OVERRIDE", "0.15"))
if (!is.finite(NOISE) || NOISE <= 0) NOISE <- 0.15

results_dir <- file.path(root_path, "experiments", "halfmoon", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

rows <- list()
for (sd in seeds) {
  set.seed(sd)
  S <- make_halfmoon_splits(n_train = N, n_test = N, noise = NOISE, seed = sd, val_frac = 0.2)
  mods <- fit_halfmoon_models(S, seed = sd)
  ksdf <- calc_ks_summary(mods, S)
  ksdf$seed <- sd
  rows[[length(rows) + 1]] <- ksdf
}

per_seed <- do.call(rbind, rows)
if ("ks_per_dim" %in% names(per_seed)) per_seed$ks_per_dim <- NULL

# Aggregate: mean and SE over seeds of the median KS across dims
agg <- do.call(rbind, lapply(split(per_seed, per_seed$model), function(df) {
  m <- mean(df$ks_median, na.rm = TRUE)
  se <- stats::sd(df$ks_median, na.rm = TRUE) / sqrt(sum(is.finite(df$ks_median)))
  data.frame(model = df$model[1], ks_median_mean = m, ks_median_se = se, stringsAsFactors = FALSE)
}))
row.names(agg) <- NULL

utils::write.csv(per_seed, file = file.path(results_dir, "ks_halfmoon_per_seed.csv"), row.names = FALSE)
utils::write.csv(agg, file = file.path(results_dir, "ks_halfmoon_summary.csv"), row.names = FALSE)

print(agg)

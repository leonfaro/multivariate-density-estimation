args <- commandArgs(trailingOnly = TRUE)
source("00_globals.R")
grid_side_arg <- NA_integer_
no_cache <- FALSE
cores_arg <- NA_integer_
timeout_sec_arg <- NA_integer_
abort_file <- NA_character_
order_arg <- NA_character_
for (a in args) {
  if (grepl("^--grid_side=", a)) grid_side_arg <- as.integer(sub("^--grid_side=", "", a))
  if (a == "--no_cache") no_cache <- TRUE
  if (grepl("^--cores=", a)) cores_arg <- as.integer(sub("^--cores=", "", a))
  if (grepl("^--timeout_sec=", a)) timeout_sec_arg <- as.integer(sub("^--timeout_sec=", "", a))
  if (grepl("^--abort_file=", a)) abort_file <- sub("^--abort_file=", "", a)
  if (grepl("^--order=", a)) order_arg <- sub("^--order=", "", a)
}
SEED <- 7L
N <- 100L
# Quick edit: set default noise here
NOISE <- 0.15
# Optional: allow override via environment variable NOISE_OVERRIDE
nv <- suppressWarnings(as.numeric(Sys.getenv("NOISE_OVERRIDE", "")))
if (is.finite(nv) && !is.na(nv) && nv > 0) NOISE <- nv
Sys.setenv(
  DATASET = "halfmoon2d",
  SEED = as.character(SEED),
  N_TRAIN = as.character(N),
  N_TEST = as.character(N),
  NOISE = as.character(NOISE)
)
source("scripts/halfmoon_data.R")
source("scripts/halfmoon_plot.R")
source("04_evaluation.R")
source("models/true_model.R")
source("models/trtf_model.R")
# Create splits and persist for reproducibility
S <- make_halfmoon_splits(n_train = N, n_test = N, noise = NOISE, seed = SEED, val_frac = 0.2)
dir.create("results", showWarnings = FALSE)
saveRDS(S, sprintf("results/splits_halfmoon2d_seed%03d.rds", SEED))
# Fit models and write NLL CSV via evaluation helper (no plotting)
mods <- fit_halfmoon_models(S, seed = SEED, order_mode = if (!is.na(order_arg)) order_arg else "as-is")
tab <- eval_halfmoon(mods, S)
# Format: round to 2 decimals; include ± 2*SE for joint
fmt2 <- function(x) sprintf("%.2f", round(as.numeric(x), 2))
if ("mean_joint_nll" %in% names(tab) && "se_joint" %in% names(tab)) {
  tab$mean_joint_nll <- sprintf("%s ± %s", fmt2(tab$mean_joint_nll), fmt2(2 * tab$se_joint))
  tab$se_joint <- NULL
}
for (nm in names(tab)) {
  if (grepl("^per_dim_nll_\\d+$", nm)) tab[[nm]] <- fmt2(tab[[nm]])
}
if ("train_test_policy" %in% names(tab)) tab$train_test_policy <- NULL
csv <- sprintf("results/nll_halfmoon_seed%03d.csv", SEED)
cat(sprintf("n=%d (Half-moon)\n", N))
print(tab)
message("Half-moon NLL (nats) CSV: ", csv)

# Also render the contour panels to the RStudio Plots pane and save PNG
clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))
N_tr <- nrow(S$X_tr)
grid_side <- if (!is.na(grid_side_arg)) grid_side_arg else clamp(floor(sqrt(100 * N_tr)), 80L, 200L)
default_cores <- tryCatch(min(8L, parallel::detectCores()), error = function(e) 1L)
cores <- if (!is.na(cores_arg) && cores_arg >= 1L) cores_arg else default_cores
timeout_sec <- if (!is.na(timeout_sec_arg) && timeout_sec_arg > 0L) timeout_sec_arg else NA_integer_
if (!is.na(cores)) cat(sprintf("Using %d workers\n", cores))
res <- plot_halfmoon_models(mods, S, grid_side = grid_side, save_png = TRUE,
                            show_plot = TRUE, no_cache = no_cache, cores = cores,
                            timeout_sec = timeout_sec, abort_file = abort_file,
                            order_mode = if (!is.na(order_arg)) order_arg else "as-is")
if (!is.null(res$png)) message("Panel PNG: ", res$png)

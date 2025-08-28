args <- commandArgs(trailingOnly = TRUE)
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
SEED <- 7L; N <- 100L; NOISE <- 0.15
Sys.setenv(
  DATASET = "halfmoon2d",
  SEED = as.character(SEED),
  N_TRAIN = as.character(N),
  N_TEST = as.character(N),
  NOISE = as.character(NOISE)
)
source("main.R"); tab <- main()
csv <- sprintf("results/nll_halfmoon_seed%03d.csv", SEED)
stopifnot(file.exists(csv))
S <- readRDS(sprintf("results/splits_halfmoon2d_seed%03d.rds", SEED))
clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))
N_tr <- nrow(S$X_tr)
grid_side <- if (!is.na(grid_side_arg)) grid_side_arg else clamp(floor(sqrt(100 * N_tr)), 80L, 200L)
default_cores <- tryCatch(min(10L, parallel::detectCores()), error = function(e) 1L)
cores <- if (!is.na(cores_arg) && cores_arg >= 1L) cores_arg else default_cores
timeout_sec <- if (!is.na(timeout_sec_arg) && timeout_sec_arg > 0L) timeout_sec_arg else NA_integer_
if (!is.na(cores)) cat(sprintf("Using %d workers\n", cores))
source("scripts/halfmoon_plot.R"); mods <- fit_halfmoon_models(S, seed = SEED,
                                                                order_mode = if (!is.na(order_arg)) order_arg else "as-is")
res <- plot_halfmoon_models(mods, S, grid_side = grid_side, save_png = TRUE,
                            show_plot = FALSE, no_cache = no_cache, cores = cores,
                            timeout_sec = timeout_sec, abort_file = abort_file,
                            order_mode = if (!is.na(order_arg)) order_arg else "as-is")
png_file <- res$png
stopifnot(file.exists(png_file))
stopifnot(identical(tab, results_table))
print(tab)
message("Panel PNG: ", png_file)
message("Half-moon NLL (nats): ", csv)

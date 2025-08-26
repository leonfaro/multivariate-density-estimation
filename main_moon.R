args <- commandArgs(trailingOnly = TRUE)
grid_side_arg <- NA_integer_
no_cache <- FALSE
for (a in args) {
  if (grepl("^--grid_side=", a)) grid_side_arg <- as.integer(sub("^--grid_side=", "", a))
  if (a == "--no_cache") no_cache <- TRUE
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
source("scripts/halfmoon_plot.R"); mods <- fit_halfmoon_models(S, seed = SEED)
res <- plot_halfmoon_models(mods, S, grid_side = grid_side, save_png = TRUE,
                            show_plot = FALSE, no_cache = no_cache)
png_file <- res$png
stopifnot(file.exists(png_file))
stopifnot(identical(tab, results_table))
print(tab)
message("Panel PNG: ", png_file)
message("Half-moon NLL (nats): ", csv)

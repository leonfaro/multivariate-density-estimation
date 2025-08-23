seed  <- 7L
ntr   <- 80L
nte   <- 80L
noise <- 0.15
Sys.setenv(
  DATASET = "halfmoon2d",
  N_TRAIN = as.character(ntr),
  N_TEST  = as.character(nte),
  NOISE   = as.character(noise),
  SEED    = as.character(seed)
)
source("main.R")
tab <- main()
csv <- sprintf("results/nll_halfmoon_seed%03d.csv", seed)
stopifnot(file.exists(csv))
S <- readRDS(sprintf("results/splits_halfmoon2d_seed%03d.rds", seed))
source("scripts/halfmoon_plot.R")
mods <- fit_halfmoon_models(S, seed = seed)
plot_halfmoon_models(mods, S, grid_n = 120, save_png = TRUE)
png_file <- sprintf("results/halfmoon_panels_seed%03d.png", seed)
stopifnot(file.exists(png_file))
stopifnot(identical(tab, results_table))
print(tab)
message("Panel PNG: ", png_file)
message("Half-moon NLL (nats): ", csv)

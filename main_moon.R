SEED <- 7L; N <- 80L; NOISE <- 0.15
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
source("scripts/halfmoon_plot.R"); mods <- fit_halfmoon_models(S, seed = SEED)
plot_halfmoon_models(mods, S, save_png = TRUE, show_plot = FALSE)
png_file <- sprintf("results/halfmoon_panels_seed%03d.png", SEED)
stopifnot(file.exists(png_file))
stopifnot(identical(tab, results_table))
print(tab)

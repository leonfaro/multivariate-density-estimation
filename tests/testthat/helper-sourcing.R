# Helper to source project code for tests, independent of working dir

root_path <- getwd()
if (basename(root_path) == "testthat") {
  root_path <- dirname(dirname(root_path))
}

source(file.path(root_path, "00_globals.R"))
source(file.path(root_path, "01_data_generation.R"))
source(file.path(root_path, "02_split.R"))
source(file.path(root_path, "R/ttm_bases.R"))
source(file.path(root_path, "R/ttm_core.R"))
source(file.path(root_path, "R/ttm_marginal.R"))
source(file.path(root_path, "R/ttm_separable.R"))
source(file.path(root_path, "R/ttm_crossterm.R"))
source(file.path(root_path, "models/true_model.R"))
source(file.path(root_path, "models/true_joint_model.R"))

# Access the Config-4D definition without triggering main()
source(file.path(root_path, "main.R"))

stderr <- function(x) stats::sd(x)/sqrt(length(x))

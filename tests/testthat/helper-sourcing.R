# Helper to source project code for tests, independent of working dir

# Robustly locate repo root by climbing up until we find expected files
find_root <- function(start = getwd()) {
  p <- normalizePath(start, winslash = "/", mustWork = FALSE)
  # In testthat, getwd() can be a temp dir; climb upwards to find the repo
  for (i in 1:10) {
    if (file.exists(file.path(p, "R", "ttm_bases.R")) &&
        file.exists(file.path(p, "ALGORITHM_SPEC.md"))) {
      return(p)
    }
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  # Fallback to testthat::test_path("..") if available
  if (exists("test_path", where = asNamespace("testthat"), inherits = FALSE)) {
    p2 <- tryCatch(testthat::test_path(".."), error = function(e) NULL)
    if (!is.null(p2)) return(normalizePath(p2, winslash = "/", mustWork = FALSE))
  }
  # Last resort: original working directory
  normalizePath(start, winslash = "/", mustWork = FALSE)
}

root_path <- find_root()

source(file.path(root_path, "00_globals.R"))
source(file.path(root_path, "01_data_generation.R"))
source(file.path(root_path, "02_split.R"))
source(file.path(root_path, "models/ttm/ttm_bases.R"))
source(file.path(root_path, "models/ttm/ttm_core.R"))
source(file.path(root_path, "models/ttm/ttm_marginal.R"))
source(file.path(root_path, "models/ttm/ttm_separable.R"))
source(file.path(root_path, "models/ttm/ttm_crossterm.R"))
source(file.path(root_path, "models/true_model.R"))
source(file.path(root_path, "models/true_joint_model.R"))

# Access the Config-4D definition without triggering main()
source(file.path(root_path, "main.R"))

stderr <- function(x) stats::sd(x)/sqrt(length(x))

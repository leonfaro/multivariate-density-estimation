# Repository loader and path helpers

repo_root <- function(start = getwd()) {
  p <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (file.exists(file.path(p, "R", "globals.R"))) return(p)
    if (file.exists(file.path(p, "code", "R", "globals.R"))) return(normalizePath(file.path(p, "code"), winslash = "/", mustWork = TRUE))
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  normalizePath(start, winslash = "/", mustWork = FALSE)
}

source_repo_modules <- function(root = repo_root()) {
  src <- function(...) source(file.path(root, ...), chdir = FALSE)
  src("R", "globals.R")
  src("R", "data_generation.R")
  src("R", "split.R")
  src("R", "evaluation.R")
  src("R", "models", "true_model.R")
  src("R", "models", "true_joint_model.R")
  src("R", "models", "trtf_model.R")
  src("R", "models", "ttm_core.R")
  src("R", "models", "ttm_marginal.R")
  src("R", "models", "ttm_separable.R")
  if (file.exists(file.path(root, "R", "models", "copula_np.R"))) {
    src("R", "models", "copula_np.R")
  }
  invisible(TRUE)
}

# Convenience for scripts that expect `root_path`
initialize_repo <- function() {
  root <- repo_root()
  assign("root_path", root, envir = .GlobalEnv)
  source_repo_modules(root)
  invisible(root)
}

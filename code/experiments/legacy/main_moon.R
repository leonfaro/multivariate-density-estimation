#!/usr/bin/env Rscript

# Half-moon runner without artifacts/logs. Uses existing utilities but disables saving.
Sys.setenv(NO_ARTIFACTS = "1")
options(mde.no_artifacts = TRUE)

args <- commandArgs(trailingOnly = TRUE)
grid_side_arg <- NA_integer_
cores_arg <- NA_integer_
timeout_sec_arg <- NA_integer_
order_arg <- NA_character_
for (a in args) {
  if (grepl("^--grid_side=", a)) grid_side_arg <- as.integer(sub("^--grid_side=", "", a))
  if (grepl("^--cores=", a)) cores_arg <- as.integer(sub("^--cores=", "", a))
  if (grepl("^--timeout_sec=", a)) timeout_sec_arg <- as.integer(sub("^--timeout_sec=", "", a))
  if (grepl("^--order=", a)) order_arg <- sub("^--order=", "", a)
}

if (!exists("locate_repo_loader", inherits = TRUE)) {
  locate_repo_loader <- function() {
    detect_script_path <- function() {
      frames <- sys.frames()
      for (i in rev(seq_along(frames))) {
        fi <- frames[[i]]
        if (!is.null(fi$ofile)) {
          path <- tryCatch(normalizePath(fi$ofile, winslash = "/", mustWork = TRUE),
                          error = function(e) NA_character_)
          if (!is.na(path) && nzchar(path)) return(path)
        }
      }
      args <- commandArgs(trailingOnly = FALSE)
      file_arg <- args[grepl("^--file=", args)]
      if (length(file_arg)) {
        cand <- sub("^--file=", "", file_arg[1])
        path <- tryCatch(normalizePath(cand, winslash = "/", mustWork = TRUE),
                         error = function(e) NA_character_)
        if (!is.na(path) && nzchar(path)) return(path)
      }
      NA_character_
    }

    start_dirs <- character()
    script_path <- detect_script_path()
    if (!is.na(script_path) && nzchar(script_path)) {
      start_dirs <- c(start_dirs, dirname(script_path))
    }
    wd <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = FALSE),
                   error = function(e) getwd())
    start_dirs <- unique(c(start_dirs, wd))
    checked <- character()
    for (start in start_dirs) {
      cur <- start
      repeat {
        cur <- tryCatch(normalizePath(cur, winslash = "/", mustWork = FALSE),
                        error = function(e) cur)
        if (!nzchar(cur) || cur %in% checked) break
        checked <- c(checked, cur)
        cand1 <- file.path(cur, "R", "loader.R")
        if (file.exists(cand1)) {
          return(normalizePath(cand1, winslash = "/", mustWork = TRUE))
        }
        cand2 <- file.path(cur, "code", "R", "loader.R")
        if (file.exists(cand2)) {
          return(normalizePath(cand2, winslash = "/", mustWork = TRUE))
        }
        parent <- dirname(cur)
        if (identical(parent, cur)) break
        cur <- parent
      }
    }
    stop("Could not locate loader.R")
  }
}

loader_path <- locate_repo_loader()
if (!exists("initialize_repo")) {
  source(loader_path, chdir = FALSE)
}
root_path <- initialize_repo()
halfmoon_dir <- file.path(root_path, "experiments", "halfmoon")
source(file.path(halfmoon_dir, "halfmoon_data.R"))
source(file.path(halfmoon_dir, "plot.R"))
source(file.path(halfmoon_dir, "true_density.R"))

# Config
SEED <- as.integer(Sys.getenv("SEED", 7))
N <- as.integer(Sys.getenv("N_MOON", 100))
NOISE <- as.numeric(Sys.getenv("NOISE", 0.15))

S <- make_halfmoon_splits(n_train = N, n_test = N, noise = NOISE, seed = SEED, val_frac = 0.2)
mods <- fit_halfmoon_models(S, seed = SEED, order_mode = if (!is.na(order_arg)) order_arg else "as-is")

# Evaluate without writing CSV
tab <- eval_halfmoon(mods, S, out_csv_path = NULL)
if ("train_test_policy" %in% names(tab)) tab$train_test_policy <- NULL
cat(sprintf("n=%d (Half-moon)\n", N))
print(tab)

# Optional contour visualization without saving PNGs
default_cores <- tryCatch(min(8L, parallel::detectCores()), error = function(e) 1L)
cores <- if (!is.na(cores_arg) && cores_arg >= 1L) cores_arg else default_cores
N_tr <- nrow(S$X_tr)
clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))
grid_side <- if (!is.na(grid_side_arg)) grid_side_arg else clamp(floor(sqrt(100 * N_tr)), 80L, 200L)
timeout_sec <- if (!is.na(timeout_sec_arg) && timeout_sec_arg > 0L) timeout_sec_arg else NA_integer_
invisible(plot_halfmoon_models(mods, S, grid_side = grid_side, save_png = FALSE,
                               show_plot = interactive(), no_cache = TRUE, cores = cores,
                               timeout_sec = timeout_sec, abort_file = NULL,
                               order_mode = if (!is.na(order_arg)) order_arg else "as-is"))

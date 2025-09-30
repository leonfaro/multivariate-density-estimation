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

args <- commandArgs(trailingOnly = TRUE)
halfmoon_dir <- file.path(root_path, "experiments", "halfmoon")
source(file.path(halfmoon_dir, "halfmoon_data.R"))
source(file.path(halfmoon_dir, "plot.R"))
source(file.path(halfmoon_dir, "true_density.R"))
grid_side_arg <- NA_integer_
no_cache <- TRUE
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
# Create splits and persist for reproducibility
S <- make_halfmoon_splits(n_train = N, n_test = N, noise = NOISE, seed = SEED, val_frac = 0.2)
# Do not persist splits to disk to avoid artifacts
# Fit models and write NLL CSV via evaluation helper (no plotting)
mods <- fit_halfmoon_models(S, seed = SEED, order_mode = if (!is.na(order_arg)) order_arg else "as-is")
timing <- attr(mods, "timing")

results_root <- file.path(halfmoon_dir, "results")
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(results_root, sprintf("seed%03d_%s", SEED, run_stamp))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

csv_path <- file.path(run_dir, sprintf("nll_halfmoon_seed%03d.csv", SEED))
eval_time <- system.time({
  tab <- eval_halfmoon(mods, S, out_csv_path = csv_path)
})[["elapsed"]]

# Record approximate prediction times per model (seconds)
timing$true$test <- system.time({ true_logdensity(S$X_te, S, Q = 32L) })[["elapsed"]]
timing$true_marg$test <- system.time({ predict(mods$true_marg, S$X_te, type = "logdensity_by_dim") })[["elapsed"]]
timing$trtf$test <- if (!is.null(mods$trtf)) {
  system.time({ predict(mods$trtf, S$X_te, type = "logdensity_by_dim") })[["elapsed"]]
} else NA_real_
timing$ttm$test <- if (!is.null(timing$ttm$test) && is.finite(timing$ttm$test)) timing$ttm$test else system.time({
  predict(mods$ttm, S$X_te, type = "logdensity_by_dim")
})[["elapsed"]]
timing$ttm_sep$test <- if (!is.null(timing$ttm_sep$test) && is.finite(timing$ttm_sep$test)) timing$ttm_sep$test else system.time({
  predict(mods$ttm_sep, S$X_te, type = "logdensity_by_dim")
})[["elapsed"]]
timing$copula_np$test <- system.time({
  predict(mods$copula_np, S$X_te, type = "logdensity_by_dim")
})[["elapsed"]]
attr(mods, "timing") <- timing

# Format: round to 2 decimals; include ± 2*SE for joint
fmt2 <- function(x) sprintf("%.2f", round(as.numeric(x), 2))
tab_fmt <- tab
if ("mean_joint_nll" %in% names(tab_fmt) && "se_joint" %in% names(tab_fmt)) {
  tab_fmt$mean_joint_nll <- sprintf("%s ± %s", fmt2(tab_fmt$mean_joint_nll), fmt2(2 * tab_fmt$se_joint))
  tab_fmt$se_joint <- NULL
}
for (nm in names(tab_fmt)) {
  if (grepl("^per_dim_nll_\\d+$", nm)) tab_fmt[[nm]] <- fmt2(tab_fmt[[nm]])
}
if ("train_test_policy" %in% names(tab_fmt)) tab_fmt$train_test_policy <- NULL
cat(sprintf("n=%d (Half-moon)\n", N))
print(tab_fmt)
message("Half-moon NLL CSV: ", csv_path)

timing$evaluation_total <- eval_time

# Also render the contour panels to the RStudio Plots pane and save PNG
clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))
N_tr <- nrow(S$X_tr)
grid_side <- if (!is.na(grid_side_arg)) grid_side_arg else clamp(floor(sqrt(3 * N_tr)), 16L, 36L)
default_cores <- tryCatch(min(8L, parallel::detectCores()), error = function(e) 1L)
cores <- if (!is.na(cores_arg) && cores_arg >= 1L) cores_arg else default_cores
timeout_sec <- if (!is.na(timeout_sec_arg) && timeout_sec_arg > 0L) timeout_sec_arg else NA_integer_
if (!is.na(cores)) cat(sprintf("Using %d workers\n", cores))
file_stub <- sprintf("halfmoon_panels_seed%03d_%s", SEED, run_stamp)
res <- plot_halfmoon_models(mods, S, grid_side = grid_side, save_png = TRUE,
                            show_plot = TRUE, no_cache = no_cache, cores = cores,
                            timeout_sec = timeout_sec, abort_file = abort_file,
                            order_mode = if (!is.na(order_arg)) order_arg else "as-is",
                            output_dir = run_dir, file_stub = file_stub)
png_path <- if (!is.null(res$png)) res$png else NULL
if (!is.null(png_path)) message("Panel PNG: ", png_path)

timing$plot_eval <- if (!is.null(res$eval_elapsed)) res$eval_elapsed else NA_real_
timing$plot_total <- if (!is.null(res$plot_elapsed)) res$plot_elapsed else NA_real_
timing$plot_png <- if (!is.null(res$png_elapsed)) res$png_elapsed else NA_real_

# Compose run summary text
summary_path <- file.path(run_dir, sprintf("halfmoon_summary_seed%03d.txt", SEED))
fmt_time <- function(x) ifelse(is.na(x) || !is.finite(x), "NA", sprintf("%.3f", x))
model_labels <- c(true = "True (Joint)", true_marg = "True (Marg)", trtf = "TRTF", ttm = "TTM Marginal",
                  ttm_sep = "TTM Separable", copula_np = "Copula NP")
timing_lines <- vapply(names(model_labels), function(nm) {
  train_t <- if (!is.null(timing[[nm]]$train)) timing[[nm]]$train else NA_real_
  test_t <- if (!is.null(timing[[nm]]$test)) timing[[nm]]$test else NA_real_
  sprintf("  %s: train=%s s, test=%s s", model_labels[[nm]], fmt_time(train_t), fmt_time(test_t))
}, character(1L), USE.NAMES = FALSE)
tab_lines <- capture.output(print(tab))
summary_lines <- c(
  sprintf("Half-moon evaluation run (%s)", run_stamp),
  sprintf("Seed: %d", SEED),
  sprintf("N_train=%d, N_val=%d, N_test=%d, noise=%.3f", S$meta$n_train, S$meta$n_val, S$meta$n_test, NOISE),
  sprintf("Grid side=%d, cores=%s", grid_side, ifelse(is.na(cores), "NA", as.character(cores))),
  sprintf("Evaluation elapsed=%.3f s", timing$evaluation_total),
  sprintf("Contour grid eval=%s s", fmt_time(timing$plot_eval)),
  sprintf("Contour plotting total=%s s, PNG write=%s s", fmt_time(timing$plot_total), fmt_time(timing$plot_png)),
  "",
  "Model timings (seconds):",
  timing_lines,
  "",
  sprintf("Evaluation table CSV: %s", csv_path),
  sprintf("Contour PNG: %s", ifelse(is.null(png_path), "(not saved)", png_path)),
  "",
  "Evaluation table (raw values):",
  tab_lines
)
writeLines(summary_lines, summary_path)
message("Summary written to ", summary_path)

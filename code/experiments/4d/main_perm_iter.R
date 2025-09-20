permutations <- list(
  c(1L, 2L, 3L, 4L),
  c(1L, 2L, 4L, 3L),
  c(1L, 3L, 2L, 4L),
  c(1L, 3L, 4L, 2L),
  c(1L, 4L, 2L, 3L),
  c(1L, 4L, 3L, 2L),
  c(2L, 1L, 3L, 4L),
  c(2L, 1L, 4L, 3L),
  c(2L, 3L, 1L, 4L),
  c(2L, 3L, 4L, 1L),
  c(2L, 4L, 1L, 3L),
  c(2L, 4L, 3L, 1L),
  c(3L, 1L, 2L, 4L),
  c(3L, 1L, 4L, 2L),
  c(3L, 2L, 1L, 4L),
  c(3L, 2L, 4L, 1L),
  c(3L, 4L, 1L, 2L),
  c(3L, 4L, 2L, 1L),
  c(4L, 1L, 2L, 3L),
  c(4L, 1L, 3L, 2L),
  c(4L, 2L, 1L, 3L),
  c(4L, 2L, 3L, 1L),
  c(4L, 3L, 1L, 2L),
  c(4L, 3L, 2L, 1L)
)

detect_script_dir <- function() {
  script_path <- tryCatch({
    normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE)
  }, error = function(e) NA_character_)
  if (is.na(script_path) || !nzchar(script_path)) {
    script_path <- tryCatch({
      args <- commandArgs(trailingOnly = FALSE)
      marker <- "--file="
      file_arg <- args[grepl(marker, args, fixed = TRUE)]
      if (length(file_arg)) {
        normalizePath(sub(marker, "", file_arg[1]), winslash = "/", mustWork = TRUE)
      } else {
        NA_character_
      }
    }, error = function(e) NA_character_)
  }
  if (!is.na(script_path) && nzchar(script_path)) {
    dirname(script_path)
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}

script_dir <- detect_script_dir()
log_file_default <- file.path(script_dir, "results", "main_perm_iter.txt")
log_file <- Sys.getenv("MAIN_PERM_ITER_LOG", log_file_default)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

seed <- as.integer(Sys.getenv("SEED", "42"))
set.seed(seed)

source(file.path(script_dir, "main_perm.R"))

start_idx_env <- suppressWarnings(as.integer(Sys.getenv("MAIN_PERM_ITER_START", "1")))
start_idx <- if (is.finite(start_idx_env) && start_idx_env >= 1L) start_idx_env else 1L
if (start_idx > length(permutations)) {
  stop(sprintf("MAIN_PERM_ITER_START=%d exceeds number of permutations (%d)",
               start_idx, length(permutations)))
}

permutations_to_run <- permutations[start_idx:length(permutations)]
total_perms <- length(permutations)
total_remaining <- length(permutations_to_run)

n_iter <- get("n", envir = .GlobalEnv)
prepare_data_original <- get("prepare_data", envir = .GlobalEnv)

set.seed(seed)
prepared_data <- prepare_data_original(n_iter, config, seed = seed)

prepare_data_cached <- function(n_arg, config_arg, seed_arg = seed) {
  prepared_data
}

assign("prepare_data", prepare_data_cached, envir = .GlobalEnv)
on.exit(assign("prepare_data", prepare_data_original, envir = .GlobalEnv), add = TRUE)

run_for_perm <- function(perm_vec, iter_index, remaining) {
  perm_string <- paste(perm_vec, collapse = ",")
  message(sprintf("[ITER] Running permutation %s (iteration %d/%d, %d remaining)",
                  perm_string, iter_index, total_perms, remaining))

  old_perm_override <- tryCatch(get("perm_override", envir = .GlobalEnv), error = function(e) NULL)
  assign("perm_override", perm_vec, envir = .GlobalEnv)

  sink(log_file, append = TRUE)
  on.exit({
    sink()
    if (!is.null(old_perm_override)) {
      assign("perm_override", old_perm_override, envir = .GlobalEnv)
    }
  }, add = TRUE)
  cat(sprintf("\n=== Permutation %s ===\n", perm_string))

  invisible(run_main_perm())
}

for (i in seq_along(permutations_to_run)) {
  perm <- permutations_to_run[[i]]
  iter_index <- start_idx - 1L + i
  remaining <- total_remaining - i + 1L
  tryCatch(run_for_perm(perm, iter_index, remaining), error = function(e) {
    message(sprintf("[ITER] Permutation %s failed: %s", paste(perm, collapse = ","), e$message))
    cat(sprintf("\n[ERROR] Permutation %s: %s\n", paste(perm, collapse = ","), e$message),
        file = log_file, append = TRUE)
  })
}

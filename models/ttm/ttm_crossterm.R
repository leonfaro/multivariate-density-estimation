## Bridge loader for Cross-term TTM implementation (repo-local, no proxies)

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

if (!exists("fit_ttm_crossterm")) {
  # Try multiple candidate paths to be robust against varying working dirs
  cand <- character(0)
  # 1) Relative to current working directory
  cand <- c(cand, file.path("models", "ttm", "fit_ttm_crossterm.R"))
  # 2) If this file is sourced by absolute path, use its directory
  this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
  if (is.character(this_file) && nzchar(this_file) && !is.na(this_file)) {
    cand <- c(cand, file.path(dirname(this_file), "fit_ttm_crossterm.R"))
  }
  # 3) If a root_path variable is defined (used elsewhere in repo)
  if (exists("root_path")) {
    cand <- c(cand, file.path(root_path, "models", "ttm", "fit_ttm_crossterm.R"))
  }
  # 4) Walk up a few parents to find repo root containing 00_globals.R
  p <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:8) {
    if (file.exists(file.path(p, "00_globals.R"))) {
      cand <- c(cand, file.path(p, "models", "ttm", "fit_ttm_crossterm.R"))
      break
    }
    np <- dirname(p); if (identical(np, p)) break; p <- np
  }
  # Deduplicate and source first hit
  cand <- unique(cand)
  hit <- NULL
  for (pth in cand) {
    if (file.exists(pth)) { hit <- pth; break }
  }
  if (!is.null(hit)) {
    source(hit)
    try(message(sprintf("[WIRE] sourced=%s", normalizePath(hit, winslash = "/", mustWork = FALSE))), silent = TRUE)
  } else {
    stop("fit_ttm_crossterm not found. Required file missing: models/ttm/fit_ttm_crossterm.R")
  }
}

# Provide a stable evaluator symbol that wraps the generic predict_ttm
if (!exists("predict_ttm_crossterm")) {
  predict_ttm_crossterm <- function(object, newdata, type = c("logdensity_by_dim", "logdensity")) {
    predict_ttm(object, newdata, match.arg(type))
  }
}

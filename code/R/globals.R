## Lightweight, defensive imports for optional packages
safe_library <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(do.call("library", list(pkg, character.only = TRUE)))
    TRUE
  } else {
    message(sprintf("[INFO] Optional package '%s' not installed; proceeding without it.", pkg))
    FALSE
  }
}

safe_library("dplyr")
safe_library("parallel")
safe_library("tram")
safe_library("trtf")
if (!exists("NC")) {
  # Fallback to 1 core if 'parallel' is unavailable
  detect_cores <- function() {
    if (requireNamespace("parallel", quietly = TRUE)) parallel::detectCores() else 1L
  }
  NC <- suppressWarnings(as.integer(detect_cores()))
  if (!is.finite(NC) || is.na(NC) || NC < 1L) NC <- 1L
}
options(mc.cores = NC)

softplus <- function(x) log1p(exp(x))

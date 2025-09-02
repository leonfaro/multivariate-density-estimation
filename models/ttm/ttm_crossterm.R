# Bridge loader for Cross-term TTM implementation
# This file sources the actual implementation from R/ttm_crossterm.R
# and ensures it works regardless of the current working directory.

if (!exists("fit_ttm_crossterm")) {
  # Try candidates deterministically
  cands <- character(0)
  if (exists("root_path")) {
    cands <- c(cands, file.path(root_path, "R", "ttm_crossterm.R"))
  }
  cands <- c(cands,
             file.path("R", "ttm_crossterm.R"),
             file.path("..", "R", "ttm_crossterm.R"),
             file.path("..", "..", "R", "ttm_crossterm.R"))
  path <- NULL
  for (p in cands) {
    if (file.exists(p)) { path <- p; break }
  }
  if (is.null(path)) {
    # Fallback: climb upwards to find repo root containing ALGORITHM_SPEC.md and R/ttm_crossterm.R
    p0 <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    cur <- p0
    for (i in 1:6) {
      if (file.exists(file.path(cur, "ALGORITHM_SPEC.md")) && file.exists(file.path(cur, "R", "ttm_crossterm.R"))) {
        path <- file.path(cur, "R", "ttm_crossterm.R"); break
      }
      par <- dirname(cur)
      if (identical(par, cur)) break
      cur <- par
    }
  }
  if (is.null(path) || !file.exists(path)) stop("Cannot locate R/ttm_crossterm.R from ttm_crossterm bridge.")
  source(path)
}

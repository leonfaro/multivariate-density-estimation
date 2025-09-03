## Bridge loader for Cross-term TTM implementation (repo-local, no proxies)

if (!exists("fit_ttm_crossterm")) {
  # Prefer concrete repo-local implementation
  p_local <- file.path("models", "ttm", "fit_ttm_crossterm.R")
  if (file.exists(p_local)) {
    source(p_local)
    try(message(sprintf("[WIRE] sourced=%s", normalizePath(p_local, winslash = "/", mustWork = FALSE))), silent = TRUE)
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

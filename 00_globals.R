# Global configuration for TRUE-Modell-Experimente
# Notation follows Theory.md

## Bibliotheken -------------------------------------------------------------
pkgs <- c("trtf", "kableExtra", "dplyr",
          "tibble", "testthat", "tidyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}
library(parallel)

softplus <- function(x) log1p(exp(x))


#' Initialize global parameters
#'
#' @return list with elements `n`, `config`, `seed`, `split_ratio`,
#'   `h_grid`, `model_ids` and `p_max`
setup_global <- function(cfg = config) {
  n <- 50
  seed <- 42
  split_ratio <- 0.5
  p_max <- 6
  h_grid <- seq_len(p_max)
  model_ids <- c("TRUE")

  list(
    n = n,
    config = cfg,
    seed = seed,
    split_ratio = split_ratio,
    h_grid = h_grid,
    model_ids = model_ids,
    p_max = p_max
  )
}

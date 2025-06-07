# Global configuration for TRUE-Modell-Experimente
# Notation follows Theory.md

## Bibliotheken -------------------------------------------------------------
library(parallel)
library(trtf)
library(ggplot2)
library(gridExtra)
library(kableExtra)
library(dplyr)
library(tibble)
library(testthat)
library(tidyr)

softplus <- function(x) log1p(exp(x))


#' Initialize global parameters
#'
#' @return list with elements N, config, seed, split_ratio,
#'   H_grid, model_ids and P_max
setup_global <- function() {
  N <- 500
  seed <- 42
  split_ratio <- 0.5
  P_max <- 6
  H_grid <- seq_len(P_max)
  model_ids <- c("TRUE")

  list(
    N = N,
    config = config,
    seed = seed,
    split_ratio = split_ratio,
    H_grid = H_grid,
    model_ids = model_ids,
    P_max = P_max
  )
}

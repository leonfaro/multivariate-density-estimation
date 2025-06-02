# Global configuration for TRUE-Modell-Experimente
# Notation follows Theory.md

# Default configuration für die Komponentenverteilungen
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = d$X1)),
  list(distr = "beta", parm = function(d) list(shape1 = d$X2, shape2 = 1)),
  list(distr = "gamma", parm = function(d) list(shape = d$X3, scale = 1))
)


#' Initialize global parameters
#'
#' @return list with elements N, config, seed, split_ratio,
#'   H_grid, model_ids and P_max
setup_global <- function() {
  N <- 500
  seed <- 42
  split_ratio <- 0.70
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

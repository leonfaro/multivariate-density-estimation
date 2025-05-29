setup_global <- function(config, P_max) {
  N <- 500
  seed <- 42
  split_ratio <- 0.70
  H_grid <- seq.int(1L, as.integer(P_max))
  model_ids <- c("TTM", "TRUE")
  list(
    N = N,
    config = config,
    seed = seed,
    split_ratio = split_ratio,
    H_grid = H_grid,
    model_ids = model_ids
  )
}

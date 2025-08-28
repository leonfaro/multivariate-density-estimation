#' Split data into train and test sets (Train/Test only)
#'
#' Randomly permutes rows using a reproducible seed and
#' applies a fixed 80/20 split ratio (train/test).
#' The seed must be identical across models to guarantee
#' comparable results.
#'
#' @param X numeric matrix of observations
#' @param seed integer RNG seed
#' @return list with `X_tr`, `X_te`
#' @export
SplitStruct <- function(X_tr, X_te) {
  list(X_tr = X_tr, X_te = X_te)
}

split_data <- function(X, seed) {
  stopifnot(is.matrix(X))
  N <- nrow(X)
  set.seed(seed)
  idx <- sample.int(N)
  n_tr  <- floor(0.8 * N)
  idx_tr  <- idx[seq_len(n_tr)]
  idx_te  <- idx[(n_tr + 1):N]
  SplitStruct(
    X[idx_tr , , drop = FALSE],
    X[idx_te , , drop = FALSE]
  )
}

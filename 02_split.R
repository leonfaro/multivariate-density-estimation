#' Split data into train, validation and test sets
#'
#' Randomly permutes rows using a reproducible seed and
#' applies a fixed 80/10/10 split ratio. The seed must
#' be identical across models to guarantee comparable
#' results.
#'
#' @param X numeric matrix of observations
#' @param seed integer RNG seed
#' @return list with `X_tr`, `X_val`, `X_te`
#' @export
SplitStruct <- function(X_tr, X_val, X_te) {
  list(X_tr = X_tr, X_val = X_val, X_te = X_te)
}

split_data <- function(X, seed) {
  stopifnot(is.matrix(X))
  N <- nrow(X)
  set.seed(seed)
  idx <- sample.int(N)
  n_tr  <- floor(0.8 * N)
  n_val <- floor(0.1 * N)
  idx_tr  <- idx[seq_len(n_tr)]
  idx_val <- idx[seq_len(n_val) + n_tr]
  idx_te  <- idx[(n_tr + n_val + 1):N]
  SplitStruct(
    X[idx_tr , , drop = FALSE],
    X[idx_val, , drop = FALSE],
    X[idx_te , , drop = FALSE]
  )
}

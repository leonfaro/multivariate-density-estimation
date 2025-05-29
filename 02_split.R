#' Split data into train and test sets
#'
#' This function implements Script 3 from `roadmap.md`.
#' It randomizes row order using a reproducible seed and returns
#' a list with training and test matrices.
#'
#' @param X numeric matrix of observations
#' @param split_ratio proportion of rows for training set
#' @param seed integer seed for RNG
#' @return list with elements `X_tr` and `X_te`
#' @export
train_test_split <- function(X, split_ratio, seed) {
  stopifnot(is.matrix(X), is.numeric(split_ratio), length(split_ratio) == 1)
  N <- nrow(X)
  set.seed(seed + 1L)
  idx <- sample.int(N)
  N_tr <- floor(split_ratio * N)
  X_tr <- X[idx[seq_len(N_tr)], , drop = FALSE]
  X_te <- X[idx[(N_tr + 1):N], , drop = FALSE]
  list(X_tr = X_tr, X_te = X_te)
}

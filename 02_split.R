#' Train-test split for a data matrix
#'
#' Splits `X` into training and test sets based on `split_ratio`.
#' A new random stream is used by adding one to `seed`.
#'
#' @param X Numeric matrix of observations
#' @param split_ratio Proportion of rows assigned to the training set
#' @param seed Integer random seed
#' @return list with components `X_tr` and `X_te`
#' @examples
#' G <- setup_global()
#' X <- gen_samples(G)
#' S <- train_test_split(X, 0.7, 42)
#' dim(S$X_tr)
#' dim(S$X_te)
train_test_split <- function(X, split_ratio, seed) {
  N <- nrow(X)
  set.seed(seed + 1)
  idx <- sample.int(N, N)
  N_tr <- floor(split_ratio * N)
  X_tr <- X[idx[seq_len(N_tr)], , drop = FALSE]
  X_te <- X[idx[(N_tr + 1):N], , drop = FALSE]
  list(X_tr = X_tr, X_te = X_te)
}

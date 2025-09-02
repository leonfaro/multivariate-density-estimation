# Standardization utilities (train-only)

standardize_train_only <- function(Xtr, X = NULL) {
  stopifnot(is.matrix(Xtr) || is.data.frame(Xtr))
  Xtr <- as.matrix(Xtr)
  mu <- colMeans(Xtr)
  sigma <- apply(Xtr, 2, sd) + .Machine$double.eps
  if (is.null(X)) {
    Xs <- sweep(sweep(Xtr, 2, mu, "-"), 2, sigma, "/")
  } else {
    Xs <- sweep(sweep(as.matrix(X), 2, mu, "-"), 2, sigma, "/")
  }
  list(Xs = Xs, mu = mu, sigma = sigma)
}


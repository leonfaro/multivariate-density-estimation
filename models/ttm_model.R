# Basic utilities for triangular transport maps
# Notation follows README.md and Theory.md

#' Center and scale each column of a matrix
#'
#' @param X Numeric matrix of size N x K
#' @param eps Small positive value for numerical stability
#' @return list with components `X_tilde`, `mu`, `sigma`
#' @export
standardize_data <- function(X, eps = 1e-8) {
  stopifnot(is.matrix(X))
  mu <- colMeans(X)
  sigma <- apply(X, 2, stats::sd) + eps
  x_tilde <- sweep(X, 2, mu, FUN = "-")
  x_tilde <- sweep(x_tilde, 2, sigma, FUN = "/")
  list(X_tilde = x_tilde, mu = mu, sigma = sigma)
}

#' Draw samples from the standard normal reference
#'
#' @param N Number of samples
#' @param K Dimension of each sample
#' @param seed RNG seed
#' @return Numeric matrix of dimension N x K
#' @export
sample_reference <- function(N, K, seed = 42) {
  stopifnot(length(N) == 1L, length(K) == 1L)
  set.seed(seed)
  matrix(rnorm(N * K), nrow = N, ncol = K)
}

#' Permute ordering of dimensions
#'
#' @param K Dimension of the space
#' @param seed RNG seed
#' @return Integer vector giving a permutation of 1:K
#' @export
shuffle_ordering <- function(K, seed = 42) {
  stopifnot(length(K) == 1L)
  set.seed(seed)
  sample.int(K, size = K, replace = FALSE)
}

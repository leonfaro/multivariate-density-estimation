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
#' @param n Number of samples
#' @param k Dimension of each sample
#' @param seed RNG seed
#' @return Numeric matrix of dimension N x K
#' @export
sample_reference <- function(n, k, seed = 42) {
  stopifnot(length(n) == 1L, length(k) == 1L)
  set.seed(seed)
  matrix(rnorm(n * k), nrow = n, ncol = k)

}

#' Permute ordering of dimensions
#'
#' @param k Dimension of the space
#' @param seed RNG seed
#' @return Integer vector giving a permutation of 1:K
#' @export
shuffle_ordering <- function(k, seed = 42) {
  stopifnot(length(k) == 1L)
  set.seed(seed)
  sample.int(k, size = k, replace = FALSE)
}

#' Container for triangular map coefficients and basis functions
#'
#' @param type Character string: "marginal", "separable" or "cross"
#' @param coeffA List of coefficient vectors for monotone terms f_k
#' @param coeffB List of coefficient vectors for offsets g_k
#' @param coeffC List of coefficient vectors for cross-terms h_k
#' @param basisF List of functions for f_k
#' @param basisG List of functions for g_k
#' @param basisH List of functions for h_k
#' @return Object of class 'MapStruct'
#' @export
MapStruct <- function(type = c("marginal", "separable", "cross"),
                      coeffA, coeffB = NULL, coeffC = NULL,
                      basisF, basisG = NULL, basisH = NULL) {
  type <- match.arg(type)
  structure(list(type = type,
                 coeffA = coeffA,
                 coeffB = coeffB,
                 coeffC = coeffC,
                 basisF = basisF,
                 basisG = basisG,
                 basisH = basisH),
            class = "MapStruct")

}

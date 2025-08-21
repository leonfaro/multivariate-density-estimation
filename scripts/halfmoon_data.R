#' Two-moons data generator and split utilities
#'
#' This script provides functions to create the classic two-moons dataset
#' and split it into training, validation and test sets.

#' Generate two-moons points
#'
#' @param n Number of samples to generate
#' @param noise Standard deviation of Gaussian noise added to the points
#' @param seed Random seed for reproducibility
#' @return Numeric matrix with columns x1 and x2
#' @export
#' @examples
#' X <- generate_two_moons(100, 0.1, 1)
#' plot(X)
generate_two_moons <- function(n, noise, seed) {
  stopifnot(is.numeric(n), n > 1, is.numeric(noise), is.numeric(seed))
  set.seed(seed)
  n1 <- floor(n / 2)
  n2 <- n - n1
  theta1 <- runif(n1, 0, pi)
  theta2 <- runif(n2, 0, pi)
  X1 <- cbind(cos(theta1), sin(theta1))
  X2 <- cbind(1 - cos(theta2), -sin(theta2) + 0.5)
  X <- rbind(X1, X2) + matrix(rnorm(2 * n, sd = noise), ncol = 2)
  idx <- sample.int(n)
  X <- X[idx, , drop = FALSE]
  colnames(X) <- c("x1", "x2")
  X
}

#' Create train/validation/test splits for two-moons data
#'
#' @param n_train Number of training samples before validation split
#' @param n_test Number of test samples
#' @param noise Noise level passed to `generate_two_moons`
#' @param seed Random seed
#' @param val_frac Fraction of training data used for validation
#' @return List with train, validation and test matrices and metadata
#' @export
make_halfmoon_splits <- function(n_train, n_test, noise, seed, val_frac = 0.2) {
  Xtr <- generate_two_moons(n_train, noise, seed)
  Xte <- generate_two_moons(n_test, noise, seed + 1)
  n_val <- max(10L, round(val_frac * n_train))
  set.seed(seed + 2)
  idx <- sample.int(n_train, n_val)
  Xval <- Xtr[idx, , drop = FALSE]
  Xtr <- Xtr[-idx, , drop = FALSE]
  S <- list(
    X_tr = Xtr,
    X_val = Xval,
    X_te = Xte,
    K = 2L,
    meta = list(
      seed = seed,
      noise = noise,
      n_train = n_train,
      n_test = n_test,
      n_val = n_val,
      val_frac = val_frac
    )
  )
  check_halfmoon_splits(S)
}

#' Validate two-moons data splits
#'
#' @param S Split list as produced by `make_halfmoon_splits`
#' @return The input list if all checks pass
#' @export
check_halfmoon_splits <- function(S) {
  stopifnot(is.list(S), all(c("X_tr", "X_val", "X_te", "K", "meta") %in% names(S)))
  mats <- c("X_tr", "X_val", "X_te")
  for (m in mats) {
    X <- S[[m]]
    if (!is.matrix(X) || ncol(X) != 2) stop(m, " must be a numeric matrix with two columns")
    if (!is.numeric(X)) stop(m, " must be numeric")
    if (any(is.na(X)) || any(is.infinite(X))) stop(m, " contains NA or Inf")
    if (is.null(colnames(X)) || any(colnames(X) != c("x1", "x2"))) stop(m, " has wrong column names")
  }
  if (S$K != 2) stop("K must be 2")
  n_train <- S$meta$n_train
  n_test <- S$meta$n_test
  n_val <- S$meta$n_val
  if (nrow(S$X_tr) + nrow(S$X_val) != n_train) stop("Training and validation sizes inconsistent")
  if (nrow(S$X_te) != n_test) stop("Test size inconsistent")
  S
}

#' Placeholder evaluation for half-moon dataset
#'
#' @param mods List of fitted models (unused)
#' @param S Split structure from `make_halfmoon_splits`
#' @param config Optional configuration (unused)
#' @return Invisible NULL after printing a short summary
#' @export
eval_halfmoon <- function(mods, S, config = NULL) {
  cat(sprintf("[HALFMOON] n_tr=%d | n_val=%d | n_te=%d\n",
              nrow(S$X_tr), nrow(S$X_val), nrow(S$X_te)))
  invisible(NULL)
}


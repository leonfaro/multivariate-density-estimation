#' Two-moons data generator and split utilities
#'
#' This script provides functions to create the classic two-moons dataset
#' and split it into training, validation and test sets.

#' Generate two-moons points
#'
#' @param n Number of samples to generate
#' @param noise Standard deviation of Gaussian noise added to the points
#' @param seed Random seed for reproducibility
#' @return List with matrix `X` (columns x1,x2) and integer vector `y` in {1,2}
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
  y <- c(rep(1L, n1), rep(2L, n2))
  X <- rbind(X1, X2) + matrix(rnorm(2 * n, sd = noise), ncol = 2)
  idx <- sample.int(n)
  X <- X[idx, , drop = FALSE]
  y <- y[idx]
  colnames(X) <- c("x1", "x2")
  list(X = X, y = y)
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
  tr <- generate_two_moons(n_train, noise, seed)
  te <- generate_two_moons(n_test, noise, seed + 1)
  Xtr <- tr$X; ytr <- tr$y
  Xte <- te$X; yte <- te$y
  n_val <- max(10L, round(val_frac * n_train))
  set.seed(seed + 2)
  idx <- sample.int(n_train, n_val)
  Xval <- Xtr[idx, , drop = FALSE]
  yval <- ytr[idx]
  Xtr <- Xtr[-idx, , drop = FALSE]
  ytr <- ytr[-idx]
  S <- list(
    X_tr = Xtr,
    X_val = Xval,
    X_te = Xte,
    y_tr = ytr,
    y_val = yval,
    y_te = yte,
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
  stopifnot(is.list(S), all(c("X_tr", "X_val", "X_te", "y_tr", "y_val", "y_te", "K", "meta") %in% names(S)))
  parts <- c("tr", "val", "te")
  for (p in parts) {
    X <- S[[paste0("X_", p)]]
    y <- S[[paste0("y_", p)]]
    if (!is.matrix(X) || ncol(X) != 2) stop(p, " X must be a numeric matrix with two columns")
    if (!is.numeric(X)) stop(p, " X must be numeric")
    if (any(is.na(X)) || any(is.infinite(X))) stop(p, " X contains NA or Inf")
    if (is.null(colnames(X)) || any(colnames(X) != c("x1", "x2"))) stop(p, " X has wrong column names")
    if (length(y) != nrow(X)) stop(p, " y length mismatch")
    if (any(is.na(y)) || any(is.infinite(y))) stop(p, " y contains NA or Inf")
    if (!all(y %in% c(1L, 2L))) stop(p, " y must contain only 1 or 2")
  }
  if (S$K != 2) stop("K must be 2")
  n_train <- S$meta$n_train
  n_test <- S$meta$n_test
  n_val <- S$meta$n_val
  if (nrow(S$X_tr) + nrow(S$X_val) != n_train) stop("Training and validation sizes inconsistent")
  if (nrow(S$X_te) != n_test) stop("Test size inconsistent")
  S
}

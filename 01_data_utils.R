# Utilities for conditional data simulation and standardisation

#' Simulate conditional sample
#'
#' @param config list of length K with entries 'distr' and 'parm'
#' @param N sample size
#' @param seed RNG seed
#' @return matrix X of dimension N x K
SimulateData <- function(config, N, seed) {
  set.seed(seed)
  K <- length(config)
  X <- matrix(NA_real_, nrow = N, ncol = K)
  for (i in seq_len(N)) {
    x_prev <- numeric(0)
    for (k in seq_len(K)) {
      ck <- config[[k]]
      pars <- EvalParams(ck$parm, x_prev)
      rfun <- get(paste0("r", ck$distr), mode = "function")
      X[i, k] <- do.call(rfun, c(list(n = 1), pars))
      x_prev <- c(x_prev, X[i, k])
    }
  }
  colnames(X) <- paste0("X", seq_len(K))
  X
}

#' Standardise matrix using train/test split
#'
#' @param X numeric matrix
#' @param rho train ratio in (0,1)
#' @return list(train = X_train, test = X_test, stats = list of (mean, sd))
Standardise <- function(X, rho) {
  N <- nrow(X)
  K <- ncol(X)
  n_train <- floor(rho * N)
  X_train <- X[seq_len(n_train), , drop = FALSE]
  X_test  <- X[seq(from = n_train + 1, to = N), , drop = FALSE]
  stats <- vector("list", K)
  for (k in seq_len(K)) {
    mu_k <- mean(X_train[, k])
    sd_k <- stats::sd(X_train[, k])
    stats[[k]] <- list(mean = mu_k, sd = sd_k)
    X_train[, k] <- (X_train[, k] - mu_k) / sd_k
    if (nrow(X_test) > 0) {
      X_test[, k] <- (X_test[, k] - mu_k) / sd_k
    }
  }
  list(train = X_train, test = X_test, stats = stats)
}

#' Evaluate parameter specification
#'
#' @param parm NULL, constant or function
#' @param z    numeric vector of previous components
#' @return evaluated parameters
EvalParams <- function(parm, z) {
  if (is.null(parm)) {
    NULL
  } else if (is.function(parm)) {
    if (length(z) == 0) {
      d <- data.frame()[, FALSE]
    } else {
      d <- as.data.frame(matrix(z, nrow = 1))
      names(d) <- paste0("X", seq_along(z))
    }
    parm(d)
  } else {
    parm
  }
}


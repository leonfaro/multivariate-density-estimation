# Triangular Transport Map Model
# Notation follows Theory.md and README.md.

# polynomial feature vector for previous variables
.poly_feat <- function(x_prev, h) {
  if (length(x_prev) == 0L) {
    return(1)
  }
  c(1, unlist(lapply(x_prev, function(x) x^(seq_len(h)))))
}

# total number of parameters for given degree and dimension
.theta_length <- function(h, K) {
  sum(2 * (1 + (0:(K - 1)) * h))
}

# unpack flat parameter vector into list of blocks
.unpack_theta <- function(theta, h, K) {
  idx <- 1L
  res <- vector("list", K)
  for (k in seq_len(K)) {
    p_len <- 1 + (k - 1) * h
    g_coef <- theta[idx:(idx + p_len - 1)]
    idx <- idx + p_len
    A_coef <- theta[idx:(idx + p_len - 1)]
    idx <- idx + p_len
    res[[k]] <- list(g = g_coef, A = A_coef)
  }
  res
}

# single-row map
.S_single <- function(x, theta_list, h) {
  K <- length(theta_list)
  z <- numeric(K)
  for (k in seq_len(K)) {
    x_prev <- if (k > 1) x[seq_len(k - 1)] else numeric(0)
    feat <- .poly_feat(x_prev, h)
    gk <- sum(theta_list[[k]]$g * feat)
    Ak <- sum(theta_list[[k]]$A * feat)
    int_val <- if (abs(Ak) < 1e-6) x[k] else (exp(Ak * x[k]) - 1) / Ak
    z[k] <- gk + int_val
  }
  z
}

# log-Jacobian for one observation
.logJ_single <- function(x, theta_list, h) {
  K <- length(theta_list)
  val <- 0
  for (k in seq_len(K)) {
    x_prev <- if (k > 1) x[seq_len(k - 1)] else numeric(0)
    feat <- .poly_feat(x_prev, h)
    Ak <- sum(theta_list[[k]]$A * feat)
    val <- val + Ak * x[k]
  }
  val
}

# negative log-likelihood for parameter vector
.neg_loglik <- function(theta, X, h) {
  K <- ncol(X)
  theta_list <- .unpack_theta(theta, h, K)
  ll <- apply(X, 1, function(x) {
    z <- .S_single(x, theta_list, h)
    log_phi <- sum(dnorm(z, log = TRUE))
    log_jac <- .logJ_single(x, theta_list, h)
    log_phi + log_jac
  })
  -mean(ll)
}

#' Fit Triangular Transport Map
#'
#' Implements Script 4 from `roadmap.md`.
#'
#' @param X_tr training matrix
#' @param X_te test matrix
#' @param H_grid integer vector of polynomial degrees
#' @return list with components `theta`, `h`, `logL_te`
#' @export
fit_TTM <- function(X_tr, X_te, H_grid) {
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  K <- ncol(X_tr)
  logL_te <- numeric(length(H_grid))
  theta_col <- vector("list", length(H_grid))

  for (i in seq_along(H_grid)) {
    h <- H_grid[i]
    init <- rep(0, .theta_length(h, K))
    opt <- optim(
      par = init,
      fn = .neg_loglik,
      X = X_tr,
      h = h,
      method = "L-BFGS-B",
      control = list(pgtol = 1e-6)
    )
    theta_col[[i]] <- opt$par
    logL_te[i] <- logL_TTM(list(theta = opt$par, h = h), X_te)
    message(sprintf("h=%s, logL_te=%f", h, logL_te[i]))
  }
  idx <- which.min(logL_te)
  list(theta = theta_col[[idx]], h = H_grid[idx], logL_te = logL_te[idx])
}

#' Negative log-likelihood for a fitted TTM model
#'
#' @param M_TTM result from `fit_TTM`
#' @param X matrix of observations
#' @return mean negative log-likelihood
#' @export
logL_TTM <- function(M_TTM, X) {
  stopifnot(is.matrix(X))
  .neg_loglik(M_TTM$theta, X, M_TTM$h)
}

# Triangular Transport Model (TTM) fitting utilities
#
# Notation follows Theory.md and README.md.

#' Negative log-likelihood of TTM for given parameters.
#'
#' @param theta numeric vector of positive slopes
#' @param X matrix of observations (rows)
#' @return mean negative log-likelihood
neg_loglik_ttm <- function(theta, X) {
  stopifnot(is.matrix(X))
  K <- ncol(X)
  Z <- sweep(X, 2, theta, FUN = "*")
  ll <- 0.5 * rowSums(Z^2) - sum(log(theta)) + 0.5 * K * log(2 * pi)
  mean(ll)
}

#' Gradient of the negative log-likelihood
#'
#' @param theta numeric vector of positive slopes
#' @param X matrix of observations
#' @return gradient vector
neg_loglik_ttm_grad <- function(theta, X) {
  stopifnot(is.matrix(X))
  means_sq <- colMeans(X^2)
  theta * means_sq - 1 / theta
}

#' Fit TTM via maximum likelihood
#'
#' Implements Algorithm in `Theory.md` with hyperparameter search.
#'
#' @param X_tr training matrix
#' @param X_te test matrix
#' @param H_grid vector of hyperparameters
#' @return list with elements `theta`, `h`, `logL_te`
#' @export
fit_TTM <- function(X_tr, X_te, H_grid) {
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  logL_te <- numeric(length(H_grid))
  theta_list <- vector("list", length(H_grid))

  for (i in seq_along(H_grid)) {
    h <- H_grid[i]
    init <- rep(1, ncol(X_tr))
    opt <- optim(
      par = init,
      fn = neg_loglik_ttm,
      gr = neg_loglik_ttm_grad,
      X = X_tr,
      method = "L-BFGS-B",
      lower = rep(0, length(init))
    )
    theta_list[[i]] <- opt$par
    logL_te[i] <- logL_TTM(list(theta = opt$par), X_te)
    message(sprintf("h=%s, logL_te=%.3f", h, logL_te[i]))
  }
  idx <- which.min(logL_te)
  list(theta = theta_list[[idx]], h = H_grid[idx], logL_te = logL_te[idx])
}

#' Compute negative log-likelihood per observation
#'
#' @param M_TTM list containing element `theta`
#' @param X matrix of observations
#' @return mean negative log-likelihood
#' @export
logL_TTM <- function(M_TTM, X) {
  neg_loglik_ttm(M_TTM$theta, X)
}

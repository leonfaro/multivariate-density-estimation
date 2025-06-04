# Triangular Transport Map using linear shift and scale
# relies on helper functions from 01_data_generation.R

#' Fit linear triangular transport map
#'
#' @param X_tr training matrix
#' @param X_te test matrix
#' @return list representing the fitted model
#' @export
fit_TTM <- function(X_tr, X_te) {
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  mu <- colMeans(X_tr)
  names(mu) <- colnames(X_tr)
  Sigma <- stats::cov(X_tr)
  L <- t(chol(Sigma))
  L_inv <- solve(L)
  model <- list(mu = mu, L = L, L_inv = L_inv)
  class(model) <- "ttm"
  model$logL_te <- logL_TTM(model, X_te)
  model
}

#' Predict log-densities for new observations
#'
#' @param object ttm model
#' @param newdata matrix of observations
#' @param type "logdensity" or "logdensity_by_dim"
#' @return numeric vector or matrix of log-densities
#' @export
predict.ttm <- function(object, newdata,
                        type = c("logdensity", "logdensity_by_dim")) {
  type <- match.arg(type)
  stopifnot(is.matrix(newdata))
  Z <- t(apply(newdata, 1L, S_forward, θ = object))
  ld_norm <- dnorm(Z, log = TRUE)
  log_diag <- log(diag(object$L_inv))
  if (type == "logdensity_by_dim") {
    res <- t(t(ld_norm) + log_diag)
    return(res)
  }
  log_det <- sum(log_diag)
  rowSums(ld_norm) + log_det
}

#' Compute mean negative log-likelihood
#'
#' @param model ttm model
#' @param X matrix of observations
#' @return scalar
#' @export
logL_TTM <- function(model, X) {
  ll <- predict(model, X, type = "logdensity")
  val <- -mean(ll)
  if (!is.finite(val)) stop("log-likelihood not finite")
  val
}

#' Negative log-likelihood by dimension
#'
#' @param model ttm model
#' @param X matrix of observations
#' @return numeric vector
#' @export
logL_TTM_dim <- function(model, X) {
  ll <- predict(model, X, type = "logdensity_by_dim")
  res <- -colMeans(ll)
  if (!all(is.finite(res))) stop("log-likelihood not finite")
  res
}

#' Conditional sampling from the fitted model
#'
#' @param model ttm model
#' @param b_fix named numeric vector of fixed coordinates
#' @param m number of samples
#' @return matrix of samples of size m x K
#' @export

#' Conditional sampling from the fitted model
#'
#' @param model ttm model
#' @param b_fix named numeric vector of fixed coordinates
#' @param m number of samples
#' @return matrix of samples of size m x K
#' @export
sample_TTM <- function(model, b_fix, m = 1L) {
  stopifnot(is.numeric(b_fix), is.numeric(m), m > 0)
  if (is.null(names(b_fix))) {
    fix_idx <- seq_along(b_fix)
  } else {
    fix_idx <- match(names(b_fix), names(model$mu))
  }
  mu <- model$mu
  L <- model$L
  Sigma <- L %*% t(L)
  A <- setdiff(seq_along(mu), fix_idx)
  Sigma_BB <- Sigma[fix_idx, fix_idx, drop = FALSE]
  Sigma_AB <- Sigma[A, fix_idx, drop = FALSE]
  Sigma_AA <- Sigma[A, A, drop = FALSE]
  mu_B <- mu[fix_idx]
  mu_A <- mu[A]
  cond_mean <- c(mu_A + Sigma_AB %*% solve(Sigma_BB) %*% (as.numeric(b_fix) - mu_B))
  cond_cov <- Sigma_AA - Sigma_AB %*% solve(Sigma_BB) %*% t(Sigma_AB)
  cond_cov <- cond_cov + diag(1e-6, nrow(cond_cov))
  Z <- MASS::mvrnorm(n = m, mu = cond_mean, Sigma = cond_cov)
  if (m == 1) Z <- matrix(Z, nrow = 1)
  out <- matrix(NA_real_, nrow = m, ncol = length(mu))
  out[, fix_idx] <- matrix(rep(as.numeric(b_fix), each = m), nrow = m)
  out[, A] <- Z
  colnames(out) <- names(mu)
  out
}

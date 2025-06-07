# TRUE baseline model: independent MLE for each dimension
# Notation follows Theory.md and README.md.

#' Negative log-likelihood for a univariate distribution.
#'
#' @param par numeric vector of parameters
#' @param x numeric vector of observations
#' @param distr character string of distribution name
#' @return scalar negative log-likelihood
neg_loglik_uni <- function(par, x, distr) {
  if (distr == "norm") {
    mu <- par[1]; sd <- par[2]
    if (sd <= 0) return(Inf)
    -sum(dnorm(x, mean = mu, sd = sd, log = TRUE))
  } else if (distr == "exp") {
    rate <- par[1]
    if (rate <= 0) return(Inf)
    x <- pmax(x, 1e-6)
    -sum(dexp(x, rate = rate, log = TRUE))
  } else if (distr == "beta") {
    a <- par[1]; b <- par[2]
    if (a <= 0 || b <= 0) return(Inf)
    x <- pmin(pmax(x, 1e-6), 1 - 1e-6)
    -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
  } else if (distr == "gamma") {
    shape <- par[1]; scale <- par[2]
    if (shape <= 0 || scale <= 0) return(Inf)
    x <- pmax(x, 1e-6)
    -sum(dgamma(x, shape = shape, scale = scale, log = TRUE))
  } else {
    stop("Unsupported distribution")
  }
}

# helper to compute log-density vector
.log_density_vec <- function(x, distr, par) {
  if (distr == "norm") {
    dnorm(x, mean = par[1], sd = par[2], log = TRUE)
  } else if (distr == "exp") {
    x <- pmax(x, 1e-6)
    dexp(x, rate = par[1], log = TRUE)
  } else if (distr == "beta") {
    x <- pmin(pmax(x, 1e-6), 1 - 1e-6)
    dbeta(x, shape1 = par[1], shape2 = par[2], log = TRUE)
  } else if (distr == "gamma") {
    x <- pmax(x, 1e-6)
    dgamma(x, shape = par[1], scale = par[2], log = TRUE)
  } else {
    stop("Unsupported distribution")
  }
}

# starting values based on method of moments
.start_par <- function(x, distr) {
  if (distr == "norm") {
    c(mean(x), sd(x))
  } else if (distr == "exp") {
    rate <- 1 / mean(x)
    if (!is.finite(rate) || rate <= 0) rate <- 1
    c(rate)
  } else if (distr == "beta") {
    m <- mean(x); v <- var(x)
    common <- m * (1 - m) / v - 1
    a <- m * common; b <- (1 - m) * common
    if (!is.finite(a) || a <= 0) a <- 1
    if (!is.finite(b) || b <= 0) b <- 1
    c(a, b)
  } else if (distr == "gamma") {
    m <- mean(x); v <- var(x)
    shape <- m^2 / v; scale <- v / m
    if (!is.finite(shape) || shape <= 0) shape <- 1
    if (!is.finite(scale) || scale <= 0) scale <- 1
    c(shape, scale)
  } else {
    stop("Unsupported distribution")
  }
}

#' Fit TRUE model via univariate MLE
#'
#' @param X_tr training matrix
#' @param X_te test matrix
#' @param config configuration list as in `setup_global()`
#' @return list with elements `theta` (list) and `logL_te`
#' @export
fit_TRUE <- function(X_tr, X_te, config) {
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  K <- length(config)
  theta_list <- vector("list", K)
  for (k in seq_len(K)) {
    distr_k <- config[[k]]$distr
    x_k <- X_tr[, k]
    init <- .start_par(x_k, distr_k)
    opt <- optim(
      par = init,
      fn = neg_loglik_uni,
      x = x_k,
      distr = distr_k,
      method = "L-BFGS-B",
      lower = rep(1e-6, length(init))
    )
    theta_list[[k]] <- opt$par
  }
  model <- list(theta = theta_list, config = config)
  logL_te <- logL_TRUE(model, X_te)
  model$logL_te <- logL_te
  model
}

#' Compute negative log-likelihood per observation
#'
#' @param M_TRUE list with element `theta`
#' @param X matrix of observations
#' @return mean negative log-likelihood
#' @export
logL_TRUE <- function(M_TRUE, X) {
  stopifnot(is.matrix(X))
  theta_list <- M_TRUE$theta
  config <- M_TRUE$config
  K <- length(theta_list)
  ll <- matrix(0, nrow = nrow(X), ncol = K)
  for (k in seq_len(K)) {
    distr_k <- config[[k]]$distr
    ll[, k] <- .log_density_vec(X[, k], distr_k, theta_list[[k]])
  }
  val <- -mean(rowSums(ll))
  if (!is.finite(val)) stop("log-likelihood not finite")
  val
}

#' Negative Log-Likelihood pro Dimension
#'
#' Berechnet f\u00fcr jedes Merkmal die durchschnittliche negative
#' Log-Likelihood unter dem TRUE-Modell.
#'
#' @param M_TRUE Liste aus `fit_TRUE`
#' @param X Matrix von Beobachtungen
#' @return numerischer Vektor mit L\u00e4nge `K`
#' @export
logL_TRUE_dim <- function(M_TRUE, X) {
  stopifnot(is.matrix(X))
  theta_list <- M_TRUE$theta
  config <- M_TRUE$config
  K <- length(theta_list)
  res <- numeric(K)
  for (k in seq_len(K)) {
    distr_k <- config[[k]]$distr
    ll_k <- .log_density_vec(X[, k], distr_k, theta_list[[k]])
    val <- -mean(ll_k)
    if (!is.finite(val)) stop("log-likelihood not finite")
    res[k] <- val
  }
  res
}

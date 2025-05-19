## Setup basic numerical utilities adhering to Notation.md --------------------
suppressPackageStartupMessages({
  library(extraDistr)
  library(tram)
  library(trtf)
})

## clip ----------------------------------------------------------------------
## Clamp values to a finite interval to avoid instabilities.
clip <- function(x, lo = 1e-6, hi = 1e6) {
  stopifnot(is.numeric(x), is.finite(lo), is.finite(hi))
  res <- pmin(hi, pmax(lo, x))
  stopifnot(all(is.finite(res)))
  res
}

## Machine epsilon (approx. 1e-10) for log-stability
EPS <- 1e-10

## numerically stable log-sum-exp ------------------------------------------------
## Compute \( \log\sum_i e^{x_i} \) without overflow.
logsumexp <- function(x) {
  stopifnot(is.numeric(x))
  m <- max(x)
  res <- m + log(sum(exp(x - m)))
  stopifnot(is.finite(res))
  res
}

## Numerically safe log of a pdf value
safe_logpdf <- function(val, eps = EPS) {
  res <- log(pmax(val, eps))
  stopifnot(all(is.finite(res)))
  res
}

## Numerically safe log-density after clipping
safe_logdens <- function(val, eps = EPS, hi = 1e6) {
  res <- log(pmin(pmax(val, eps), hi))
  stopifnot(all(is.finite(res)))
  res
}

## log-sum-exp for stable aggregation
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

## Bound a CDF away from 0 and 1
safe_cdf <- function(val, eps = EPS) {
  res <- pmax(eps, pmin(1 - eps, val))
  stopifnot(all(is.finite(res)))
  res
}

## Default configuration for three illustrative distributions
config3 <- list(
  list(distr = "norm",  parm = NULL),
  list(distr = "exp",   parm = function(d) list(rate = exp(d$X1))),
  list(
    distr = "gamma",
    parm  = function(d) list(
      shape = d$X2^2,
      rate  = 1)
    )
  )
config <- config3

## dimension K of the target random vector X_pi
K <- length(config)

SEED <- 24
suppressPackageStartupMessages({
  library(extraDistr)
  library(tram)
  library(trtf)
})

softplus <- function(x) {
  stopifnot(is.numeric(x))
  ifelse(x > 20, x, log1p(exp(x)))
}

EPS <- 1e-10

safe_cdf <- function(val, eps = EPS) {
  res <- pmax(eps, pmin(1 - eps, val))
  stopifnot(all(is.finite(res)))
  res
}

config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "gamma", parm  = function(d) list(shape = softplus(d$X2), rate = softplus(d$X1)))
)


K <- length(config)

source("01_map_definition_S.R")
source("02_sampling.R")
source("03_param_baseline.R")

dist_registry <- make_dist_registry()

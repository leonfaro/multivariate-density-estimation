suppressPackageStartupMessages({
  library(extraDistr)
  library(tram)
  library(trtf)
})

clip <- function(x, lo = 1e-6, hi = 1e6) pmin(hi, pmax(lo, x))
EPS <- 1e-10
safe_logpdf <- function(val, eps = EPS) log(pmax(val, eps))
safe_cdf    <- function(val, eps = EPS) pmax(eps, pmin(1 - eps, val))

config3 <- list(
  list(distr = "norm",  parm = NULL),
  list(distr = "exp",   parm = function(d) list(rate = exp(d$X1))),
  list(distr = "gamma", parm = function(d) list(shape = d$X2, rate = 1))
)
config <- config3
K <- length(config)

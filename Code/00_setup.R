suppressPackageStartupMessages({
  library(extraDistr)
  library(tram)
  library(trtf)
  if (requireNamespace("ggplot2", quietly = TRUE))
    library(ggplot2)
})

clip <- function(x, lo = 1e-6, hi = 1e6) pmin(hi, pmax(lo, x))
EPS <- 1e-10
safe_logpdf <- function(val, eps = EPS) log(pmax(val, eps))
safe_cdf    <- function(val, eps = EPS) pmax(eps, pmin(1 - eps, val))

config4 <- list(
  list(distr = "norm",    parm = NULL),
  list(distr = "t",       parm = function(d) list(df = 3 + 0.5 * d$X1)),
  list(distr = "laplace", parm = function(d) list(m = 0.3 * d$X2,
                                                  s = clip(1 + 0.1 * d$X2))),
  list(distr = "logis",   parm = function(d) list(location = 0.2 * d$X3,
                                                  scale = clip(1 + 0.05 * d$X3)))
)
config <- config4
K <- length(config)

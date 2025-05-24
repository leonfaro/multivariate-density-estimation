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

safe_logcdf <- function(val, eps = EPS) {
  lo <- log(eps)
  hi <- log1p(-eps)
  res <- pmax(lo, pmin(hi, val))
  stopifnot(all(is.finite(res)))
  res
}

safe_pars <- function(pars, dname) {
  if (dname == "norm" && !is.null(pars$sd)) stopifnot(all(pars$sd > 0))
  if (dname == "exp"  && !is.null(pars$rate)) stopifnot(all(pars$rate > 0))
  if (dname == "gamma") {
    if (!is.null(pars$shape)) stopifnot(all(pars$shape > 0))
    if (!is.null(pars$rate))  stopifnot(all(pars$rate  > 0))
  }
  if (dname == "weibull") {
    if (!is.null(pars$shape)) stopifnot(all(pars$shape > 0))
    if (!is.null(pars$scale)) stopifnot(all(pars$scale > 0))
  }
  if (dname == "lnorm"  && !is.null(pars$sdlog)) stopifnot(all(pars$sdlog > 0))
  if (dname == "pois"   && !is.null(pars$lambda)) stopifnot(all(pars$lambda >= 0))
  if (dname %in% c("bern", "binom") && !is.null(pars$prob))
    stopifnot(all(pars$prob > 0 & pars$prob < 1))
  if (dname == "binom"  && !is.null(pars$size)) stopifnot(all(pars$size >= 1))
  if (dname == "beta") {
    if (!is.null(pars$shape1)) stopifnot(all(pars$shape1 > 0))
    if (!is.null(pars$shape2)) stopifnot(all(pars$shape2 > 0))
  }
  if (dname == "logis" && !is.null(pars$scale)) stopifnot(all(pars$scale > 0))
  pars
}

q_supports_logp <- c(
  norm    = TRUE,
  exp     = TRUE,
  gamma   = TRUE,
  weibull = TRUE,
  lnorm   = TRUE,
  pois    = TRUE,
  t       = TRUE,
  laplace = FALSE,
  beta    = TRUE,
  logis   = TRUE,
  sn      = TRUE
)

p_supports_logp <- q_supports_logp

safe_support <- function(x, dname, pars = list()) {
  valid <- switch(dname,
    exp     = x > 0,
    gamma   = x > 0,
    weibull = x > 0,
    lnorm   = x > 0,
    beta    = x > 0 & x < 1,
    pois    = x >= 0 & floor(x) == x,
    bern    = x %in% c(0, 1),
    binom   = {
      size <- if (!is.null(pars$size)) pars$size else Inf
      x >= 0 & x <= size & floor(x) == x
    },
    TRUE
  )
  if (!all(valid)) stop("value outside support for ", dname)
  x
}

if (!exists("config"))
  stop("config must be defined before sourcing 00_setup.R")

K <- length(config)

source("01_map_definition_S.R", local = TRUE)
source("02_sampling.R", local = TRUE)
source("03_baseline.R", local = TRUE)

dist_registry <- make_dist_registry()

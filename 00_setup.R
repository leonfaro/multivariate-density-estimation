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

# zugelassene Verteilungen laut Spezifikation
allowed_dists_full <- c(
  "norm", "lnorm", "t", "skewt",
  "exp", "weibull", "laplace", "gpd", "burrxii",
  "invgauss", "ncchisq", "gamma", "beta"
)

# wendet die in dist_registry definierten Linkfunktionen an
apply_links <- function(pars, dname) {
  if (dname %in% names(dist_registry)) {
    spec <- dist_registry[[dname]]
    for (j in seq_along(spec$param_names)) {
      nm <- spec$param_names[j]
      if (!is.null(pars[[nm]]))
        pars[[nm]] <- link_fns[[spec$link_vector[j]]](pars[[nm]])
    }
  }
  pars
}

# prüft erlaubte Verteilungen und Typen
validate_config <- function(cfg) {
  for (k in seq_along(cfg)) {
    ck <- cfg[[k]]
    stopifnot(is.list(ck))
    stopifnot(ck$distr %in% allowed_dists_full)
    if (!is.null(ck$parm)) {
      stopifnot(is.function(ck$parm))
      dummy <- if (k > 1) {
        setNames(as.data.frame(matrix(0, nrow = 1, ncol = k - 1)),
                 paste0("X", seq_len(k - 1)))
      } else {
        data.frame()[, FALSE]
      }
      pars <- ck$parm(dummy)
      pars <- apply_links(pars, ck$distr)
      safe_pars(pars, ck$distr)
    }
  }
  cfg
}

config <- validate_config(config)

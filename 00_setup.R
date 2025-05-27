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


link_fns <- list(identity = function(z) z, softplus = softplus)

make_dist_registry <- function() {
  reg <- list(
    norm = list(
      dim = 1,
      param_names = c("mu", "sigma"),
      link_vector = c("identity", "softplus"),
      logpdf = function(x, mu, sigma) dnorm(x, mean = mu, sd = sigma, log = TRUE),
      invcdf = qnorm
    ),
    exp = list(
      dim = 1,
      param_names = "rate",
      link_vector = "softplus",
      logpdf = function(x, rate) dexp(x, rate = rate, log = TRUE),
      invcdf = qexp
    ),
    gamma = list(
      dim = 1,
      param_names = c("shape", "rate"),
      link_vector = c("softplus", "softplus"),
      logpdf = function(x, shape, rate) dgamma(x, shape = shape, rate = rate, log = TRUE),
      invcdf = qgamma
    ),
    weibull = list(
      dim = 1,
      param_names = c("shape", "scale"),
      link_vector = c("softplus", "softplus"),
      logpdf = function(x, shape, scale) dweibull(x, shape = shape, scale = scale, log = TRUE),
      invcdf = qweibull
    ),
    lnorm = list(
      dim = 1,
      param_names = c("meanlog", "sdlog"),
      link_vector = c("identity", "softplus"),
      logpdf = function(x, meanlog, sdlog) dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE),
      invcdf = qlnorm
    ),
    pois = list(
      dim = 1,
      param_names = "lambda",
      link_vector = "softplus",
      logpdf = function(x, lambda) dpois(x, lambda = lambda, log = TRUE),
      invcdf = qpois
    ),
    beta = list(
      dim = 1,
      param_names = c("shape1", "shape2"),
      link_vector = c("softplus", "softplus"),
      logpdf = function(x, shape1, shape2) dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE),
      invcdf = qbeta
    ),
    logis = list(
      dim = 1,
      param_names = c("location", "scale"),
      link_vector = c("identity", "softplus"),
      logpdf = function(x, location, scale) dlogis(x, location = location, scale = scale, log = TRUE),
      invcdf = qlogis
    )
  )
  reg
}

if (!exists("config"))
  stop("config must be defined before sourcing 00_setup.R")

dist_registry <- make_dist_registry()

# zugelassene Verteilungen laut Spezifikation
allowed_dists_full <- c(
  "norm", "lnorm", "t", "skewt",
  "exp", "weibull", "laplace", "gpd", "burrxii",
  "invgauss", "ncchisq", "gamma", "beta", "logis"
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

dist_fun <- function(pref, name) get(paste0(pref, name))

get_pars <- function(k, x_prev, cfg) {
  ck <- cfg[[k]]
  if (is.null(ck$parm)) return(list())
  if (length(x_prev) == 0) {
    x_df <- data.frame()
  } else if (is.null(dim(x_prev))) {
    names(x_prev) <- paste0("X", seq_along(x_prev))
    x_df <- as.data.frame(as.list(x_prev))
  } else {
    colnames(x_prev) <- paste0("X", seq_len(ncol(x_prev)))
    x_df <- as.data.frame(x_prev)
  }
  pars <- ck$parm(x_df)
  pars <- lapply(pars, function(x) if (is.null(x)) 1 else x)
  pars <- apply_links(pars, ck$distr)
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
    }
  }
  cfg
}

config <- validate_config(config)
K <- length(config)

source("01_map_definition_S.R", local = TRUE)
source("02_sampling.R", local = TRUE)
source("03_baseline.R", local = TRUE)

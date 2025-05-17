# Setup and configuration for K=4 triangular transport example
# This script defines helper utilities, distribution parameters and
# numerical wrappers used throughout the K4 demos.  It should be sourced
# prior to generating samples or fitting models.
suppressPackageStartupMessages({
  library(extraDistr)  # für Laplace
  library(tram)        # für Lm, Colr, BoxCox, mlt
  library(trtf)        # optional
  library(ggplot2)     # für Plots
})

# --- Block A: Konfiguration ---
# Distributional setup for the four dimensional example.  Each list
# element contains the distribution name and a function returning the
# conditional parameters given previously sampled coordinates.
#' Clamp values into a finite interval
#'
#' Utility function used in parameter functions to keep standard deviations
#' and other scale parameters within a sensible range.
#'
#' @param x Numeric vector to clamp.
#' @param lo Lower bound.
#' @param hi Upper bound.
#' @return Numeric vector with values truncated to `[lo, hi]`.
clip <- function(x, lo = 1e-6, hi = 1e6) pmin(hi, pmax(lo, x))
EPS <- 1e-10
config4 <- list(
  list(distr = "norm",    parm = NULL),
  list(distr = "t",       parm = function(d) list(df = 3 + 0.5 * d$X1)),
  list(distr = "laplace", parm = function(d) list(m = 0.3 * d$X2, s = clip(1 + 0.1 * d$X2))),
  list(distr = "logis",   parm = function(d) list(location = 0.2 * d$X3, scale = clip(1 + 0.05 * d$X3)))
)
config <- config4; K <- length(config)

# --- Numerisch stabile Hilfsfunktionen ---
#' Helper returning distribution functions by name
#'
#' @param pref Prefix, e.g. "d" for density, "p" for cdf, "q" for quantile.
#' @param name Short distribution name ("norm", "t", ...).
#' @return Function object corresponding to the requested distribution.
dist_fun <- function(pref, name) get(paste0(pref, name))
#' Numerically safe density and cdf helpers
#'
#' These wrappers avoid taking logs of zero and keep cumulative
#' distribution values away from the boundaries.
#'
#' @param val Numeric vector of values.
#' @param eps Small positive constant used for clipping.
#' @return Numeric vector of transformed values.
safe_logpdf <- function(val, eps = EPS) log(pmax(val, eps))
safe_cdf    <- function(val, eps = EPS) pmax(eps, pmin(1 - eps, val))
#' Obtain conditional distribution parameters for dimension `k`
#'
#' @param k    Dimension index.
#' @param x_prev Previously sampled values `X1, ..., X_{k-1}`.
#' @param cfg  List describing the DGP configuration.
#' @return List of parameters for the chosen distribution.
get_pars <- function(k, x_prev, cfg) {
  ck <- cfg[[k]]
  if (is.null(ck$parm)) return(list())
  names(x_prev) <- paste0("X", seq_along(x_prev))
  ck$parm(as.data.frame(as.list(x_prev)))
}
#' Density contribution for dimension `k`
#'
#' @param k      Dimension index.
#' @param xk     Observed value for dimension `k`.
#' @param x_prev Vector of previous coordinates.
#' @param cfg    DGP configuration list.
#' @param log    Logical; return log-density if `TRUE`.
#' @return Scalar density or log-density value.
pdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  pars <- get_pars(k, x_prev, cfg)
  dens <- do.call(dist_fun("d", cfg[[k]]$distr), c(list(x = xk), pars, list(log = FALSE)))
  if (log) safe_logpdf(dens) else dens
}

#' CDF contribution for dimension `k`
#'
#' @inheritParams pdf_k
#' @return Scalar CDF or log-CDF value.
cdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  dname <- cfg[[k]]$distr
  pars  <- get_pars(k, x_prev, cfg)
  args  <- c(
    if (dname == "sn") list(x = xk) else list(q = xk),
    pars,
    list(log.p = FALSE)
  )
  cdfv <- do.call(dist_fun("p", dname), args)
  cdfv <- safe_cdf(cdfv)
  if (log) log(cdfv) else cdfv
}

#' Quantile transform step of the triangular transport
#'
#' @inheritParams pdf_k
#' @param u Univariate probability value to invert.
#' @param eps Numerical stability constant.
#' @return Numeric quantile value.
qtf_k <- function(k, u, x_prev, cfg, eps = EPS) {
  u <- pmin(1 - eps, pmax(eps, u)); dname <- cfg[[k]]$distr
  pars <- get_pars(k, x_prev, cfg)
  if (dname == "sn") {
    res <- tryCatch(
      do.call(dist_fun("q", dname), list(p = u, dp = c(unlist(pars), tau = 0), solver = "RFB")),
      error = function(e) do.call(dist_fun("q", dname), list(p = u, dp = c(unlist(pars), tau = 0), solver = "NR"))
    )
  } else {
    res <- do.call(dist_fun("q", dname), c(list(p = u), pars))
  }
  res[is.infinite(res)] <- sign(res[is.infinite(res)]) * 1e6
  res[is.na(res)]       <- 0
  res
}

# --- DGP via Triangular Transport ---
# The following helpers allow to generate data from the specified
# distributional configuration using a triangular transport map.

#' Determinant of the Jacobian of the inverse transport
#'
#' @param logd Matrix of log-density contributions for each dimension.
#' @return Numeric vector with determinants for each observation.
det_J <- function(logd) rowSums(logd)

#' Log-likelihood of a transformed sample
#'
#' @param Z_eta Latent normal variables after transport.
#' @param logdet Log-determinant returned by `det_J`.
#' @return Numeric vector of log-likelihood contributions.
loglik <- function(Z_eta, logdet) {
  -0.5 * rowSums(Z_eta^2) - (ncol(Z_eta) / 2) * log(2*pi) + logdet
}

#' Sample from the standard normal reference measure
#'
#' @param N Sample size.
#' @return List with U(0,1) and Z values.
eta_sample <- function(N) {
  Z_eta <- matrix(rnorm(N * K), nrow = N)
  U_eta <- pnorm(Z_eta)
  list(U_eta = U_eta, Z_eta = Z_eta)
}

#' Apply inverse transport S^{-1} to obtain data `X_pi`
#'
#' @param U_eta Matrix of uniform variables.
#' @param cfg   Configuration list.
#' @return List with X, U and log-density contributions.
S_inv <- function(U_eta, cfg = config) {
  U_eta <- apply(U_eta, 2, function(u) pmin(1 - EPS, pmax(EPS, u)))
  n <- nrow(U_eta); X_pi <- matrix(NA_real_, n, K); logd <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    x_prev <- numeric(0)
    for (k in seq_len(K)) x_prev[k] <- qtf_k(k, U_eta[i,k], x_prev, cfg)
    X_pi[i,] <- x_prev
    for (k in seq_len(K)) {
      logd[i,k] <- pdf_k(k, X_pi[i,k], X_pi[i,seq_len(k-1)], cfg, log = TRUE) -
                   dnorm(qnorm(U_eta[i,k]), log = TRUE)
    }
  }
  list(X_pi = X_pi, U_eta = U_eta, Z_eta = qnorm(U_eta), logd = logd)
}

#' Draw a sample from the configured multivariate distribution
#'
#' @param N   Number of observations.
#' @param cfg Configuration list.
#' @return List with samples and log-density information.
pi_sample <- function(N, cfg = config) {
  ref <- eta_sample(N); inv <- S_inv(ref$U_eta, cfg)
  list(X_pi = inv$X_pi, U_eta = ref$U_eta, Z_eta = inv$Z_eta, logd = inv$logd)
}

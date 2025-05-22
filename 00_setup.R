## Grundsetup ---------------------------------------------------------------
# kurze Helfer und Beispiel-Configs für X_pi
# Eingabe: `config_choice` (3 oder 4), Standard 3
# Ausgabe: `config`, `K`, Funktionen wie `softplus()` und `logsumexp()`
SEED <- 24
suppressPackageStartupMessages({
  library(extraDistr)
  library(tram)
  library(trtf)
})

## softplus ------------------------------------------------------------------
## Glatte Positivabbildung für Parameterwerte
softplus <- function(x) {
  stopifnot(is.numeric(x))
  ifelse(x > 50, x, log1p(exp(x)))
}

## Softplus -------------------------------------------------------------------
## strictly increasing map R->(0,Inf)
softplus <- function(x) {
  ifelse(x > 20, x, log1p(exp(x)))
}

## Maschinen-Epsilon ~1e-10 für Log-Stabilität
EPS <- 1e-10

## numerisch stabiles log-sum-exp ---------------------------------------------
## berechnet log(sum(exp(x)))
logsumexp <- function(x) {
  stopifnot(is.numeric(x))
  m <- max(x)
  res <- m + log(sum(exp(x - m)))
  stopifnot(is.finite(res))
  res
}

## numerisch sichere Log-Dichte nach Clipping
safe_logdens <- function(val, eps = EPS, hi = 1e6) {
  res <- log(pmin(pmax(val, eps), hi))
  stopifnot(all(is.finite(res)))
  res
}

## CDF weg von 0 und 1 halten
safe_cdf <- function(val, eps = EPS) {
  res <- pmax(eps, pmin(1 - eps, val))
  stopifnot(all(is.finite(res)))
  res
}

## Standard-Config mit drei Verteilungen
config3 <- list(
  ## normal with mean 0, sd 1
  list(distr = "norm", parm = NULL),
  ## exponential rate depends on previous dimension
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  ## gamma shape depends on X2, fixed rate 1
  list(distr = "gamma", parm  = function(d)
    list(shape = softplus(d$X2), rate = 1))
)

## Alternative Config mit vier Verteilungen
config4 <- list(
  list(distr = "norm",    parm = NULL),
  list(distr = "t",       parm = function(d) list(df = 3 + 0.5 * d$X1)),
  list(distr = "laplace", parm = function(d)
    list(m = 0.3 * d$X2, s = softplus(1 + 0.1 * d$X2))),
  list(distr = "logis",   parm = function(d)
    list(location = 0.2 * d$X3, scale = softplus(1 + 0.05 * d$X3)))
)

if (!exists("config_choice")) config_choice <- 3
config <- if (config_choice == 4) config4 else config3

## Dimension K des Zielvektors X_pi
K <- length(config)

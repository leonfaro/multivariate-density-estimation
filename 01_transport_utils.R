
# Hilfsfunktionen für triangulären Transport
# Eingabe: Liste `config` mit K Verteilungen
# Ausgabe: Funktionen `pdf_k`, `cdf_k`, `qtf_k`, `S_inv`, `pi_sample`
# Ablauf k=1..K
#   1. Parameter aus früheren X berechnen
#   2. U_eta durch `qtf_k` mappen, logd sammeln
# Map ist monoton, daher invertierbar

SAFE_PAR_COUNT <- 0
SAFE_SUPPORT_COUNT <- 0

clip_count <- function(x, lo, hi) {
  changed <- (x < lo) | (x > hi)
  list(val = pmin(hi, pmax(lo, x)), count = sum(changed))
}

safe_pars <- function(pars, dname) {
  lo <- EPS
  hi <- 1e6
  if (dname == "norm" && !is.null(pars$sd)) {
    res <- clip_count(pars$sd, lo, hi)
    SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
    pars$sd <- res$val
  }
  if (dname == "exp" && !is.null(pars$rate)) {
    res <- clip_count(pars$rate, lo, hi)
    SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
    pars$rate <- res$val
  }
  if (dname == "gamma") {
    if (!is.null(pars$shape)) {
      res <- clip_count(pars$shape, lo, hi)
      SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
      pars$shape <- res$val
    }
    if (!is.null(pars$rate)) {
      res <- clip_count(pars$rate, lo, hi)
      SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
      pars$rate <- res$val
    }
  }

  if (dname == "weibull") {
    if (!is.null(pars$shape)) {
      res <- clip_count(pars$shape, lo, hi)
      SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
      pars$shape <- res$val
    }
    if (!is.null(pars$scale)) {
      res <- clip_count(pars$scale, lo, hi)
      SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
      pars$scale <- res$val
    }
  }
  if (dname == "lnorm" && !is.null(pars$sdlog)) {
    res <- clip_count(pars$sdlog, lo, hi)
    SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
    pars$sdlog <- res$val
  }
  if (dname == "pois" && !is.null(pars$lambda)) {
    res <- clip_count(pars$lambda, lo, hi)
    SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
    pars$lambda <- res$val
  }
  if (dname %in% c("bern", "binom") && !is.null(pars$prob)) {
    changed <- (pars$prob < lo) | (pars$prob > 1 - lo)
    SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + sum(changed)
    pars$prob <- pmin(1 - lo, pmax(lo, pars$prob))
  }
  if (dname == "binom" && !is.null(pars$size)) {
    res <- clip_count(pars$size, 1, hi)
    SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
    pars$size <- res$val
  }
  if (dname == "beta") {
    if (!is.null(pars$shape1)) {
      res <- clip_count(pars$shape1, lo, hi)
      SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
      pars$shape1 <- res$val
    }
    if (!is.null(pars$shape2)) {
      res <- clip_count(pars$shape2, lo, hi)
      SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
      pars$shape2 <- res$val
    }
  }
  if (dname == "logis" && !is.null(pars$scale)) {
    res <- clip_count(pars$scale, lo, hi)
    SAFE_PAR_COUNT <<- SAFE_PAR_COUNT + res$count
    pars$scale <- res$val
  }
  pars
}

safe_support <- function(x, dname, pars = list()) {
  lo <- EPS
  hi <- 1e6
  old <- x
  res <- switch(dname,
    exp     = clip(x, lo, hi),
    gamma   = clip(x, lo, hi),
    weibull = clip(x, lo, hi),
    lnorm   = clip(x, lo, hi),
    beta    = pmin(1 - lo, pmax(lo, x)),
    pois    = pmax(0, round(x)),
    bern    = ifelse(x > 0.5, 1, 0),
    binom   = {
      size <- if (!is.null(pars$size)) pars$size else hi
      pmin(size, pmax(0, round(x)))
    },
    x
  )
  SAFE_SUPPORT_COUNT <<- SAFE_SUPPORT_COUNT + sum(old != res)
  res
}


dist_fun <- function(pref, name) get(paste0(pref, name))
get_pars <- function(k, x_prev, cfg) {
  ck <- cfg[[k]]
  if (is.null(ck$parm)) return(list())
  if (is.null(dim(x_prev))) {
    names(x_prev) <- paste0("X", seq_along(x_prev))
    x_df <- as.data.frame(as.list(x_prev))
  } else {
    colnames(x_prev) <- paste0("X", seq_len(ncol(x_prev)))
    x_df <- as.data.frame(x_prev)
  }
  pars <- ck$parm(x_df)
  safe_pars(pars, ck$distr)
}

pdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  pars <- get_pars(k, x_prev, cfg)
  xk   <- safe_support(xk, cfg[[k]]$distr, pars)
  dens <- do.call(dist_fun("d", cfg[[k]]$distr),
                  c(list(x = xk), pars, list(log = log)))
  dens
}

cdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  dname <- cfg[[k]]$distr
  pars  <- get_pars(k, x_prev, cfg)
  xk    <- safe_support(xk, dname, pars)
  args  <- c(if (dname == "sn") list(x = xk) else list(q = xk),
             pars,
             list(log.p = FALSE))
  cdfv <- do.call(dist_fun("p", dname), args)
  cdfv <- safe_cdf(cdfv)
  if (log) log(cdfv) else cdfv
}

qtf_k <- function(k, u, x_prev, cfg, log.p = FALSE) {
  dname <- cfg[[k]]$distr
  pars <- get_pars(k, x_prev, cfg)
  if (dname == "sn") {
    res <- tryCatch(
      do.call(dist_fun("q", dname),
              list(p = u, dp = c(unlist(pars), tau = 0), solver = "RFB",
                   log.p = log.p)),
      error = function(e)
        do.call(dist_fun("q", dname),
                list(p = u, dp = c(unlist(pars), tau = 0), solver = "NR",
                     log.p = log.p))
    )
  } else {
    res <- do.call(dist_fun("q", dname), c(list(p = u, log.p = log.p), pars))
  }
  safe_support(res, dname, pars)
}

logdet_J <- function(logd) rowSums(logd)
loglik <- function(Z_eta, logdet_J) {
  -0.5 * rowSums(Z_eta^2) - (ncol(Z_eta)/2) * log(2 * pi) + logdet_J
}
eta_sample <- function(N) {
  Z_eta <- matrix(rnorm(N * K), nrow = N)
  U_eta <- pnorm(Z_eta)
  list(U_eta = U_eta, Z_eta = Z_eta)
}
S_inv <- function(U_eta, cfg = config, Z_eta = qnorm(U_eta)) {
  n <- nrow(U_eta)
  X_pi <- matrix(NA_real_, n, K)
  logd <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    x_prev <- numeric(0)
    for (k in seq_len(K)) {
      logu <- pnorm(Z_eta[i, k], log.p = TRUE)
      x_prev[k] <- qtf_k(k, logu, x_prev, cfg, log.p = TRUE)
    }
    X_pi[i, ] <- x_prev
    for (k in seq_len(K)) {
      logd[i, k] <- pdf_k(k, X_pi[i, k], X_pi[i, seq_len(k-1)], cfg, log = TRUE) -
                    dnorm(Z_eta[i, k], log = TRUE)
    }
  }
  list(X_pi = X_pi, U_eta = U_eta, Z_eta = Z_eta, logd = logd)
}
pi_sample <- function(N, cfg = config) {
  ref <- eta_sample(N)
  inv <- S_inv(ref$U_eta, cfg, ref$Z_eta)
  list(X_pi = inv$X_pi, U_eta = ref$U_eta,
       Z_eta = inv$Z_eta, logd = inv$logd)
}

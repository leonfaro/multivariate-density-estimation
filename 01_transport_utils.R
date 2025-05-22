
# Hilfsfunktionen für triangulären Transport
# Eingabe: Liste `config` mit K Verteilungen
# Ausgabe: Funktionen `pdf_k`, `cdf_k`, `qtf_k`, `S_inv`, `pi_sample`
# Ablauf k=1..K
#   1. Parameter aus früheren X berechnen
#   2. U_eta durch `qtf_k` mappen, logd sammeln
# Map ist monoton, daher invertierbar


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

## capability table for log.p in quantiles
q_supports_logp <- c(
  norm    = TRUE,
  exp     = TRUE,
  gamma   = TRUE,
  pois    = TRUE,
  t       = TRUE,
  laplace = FALSE,
  logis   = TRUE,
  sn      = TRUE
)

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
  support <- q_supports_logp[dname]
  if (is.na(support))
    stop("log.p capability not specified for distribution ", dname)
  if (log.p && !support) {
    u <- exp(u)
    log.p <- FALSE
  }
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

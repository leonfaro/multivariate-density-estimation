
# Hilfsfunktionen für triangulären Transport
# Eingabe: Liste `config` mit K Verteilungen
# Ausgabe: Funktionen `pdf_k`, `cdf_k`, `qtf_k`, `S_inv`, `pi_sample`
# Ablauf k=1..K
#   1. Parameter aus früheren X berechnen
#   2. U_eta durch `qtf_k` mappen, logd sammeln
# Map ist monoton, daher invertierbar



safe_pars <- function(pars, dname) {
  if (dname == "norm" && !is.null(pars$sd)) {
    pars$sd <- softplus(pars$sd)
  }
  if (dname == "exp" && !is.null(pars$rate)) {
    pars$rate <- softplus(pars$rate)
  }
  if (dname == "gamma") {
    if (!is.null(pars$shape)) pars$shape <- softplus(pars$shape)
    if (!is.null(pars$rate)) pars$rate <- softplus(pars$rate)
  }
  if (dname == "weibull") {
    if (!is.null(pars$shape)) pars$shape <- softplus(pars$shape)
    if (!is.null(pars$scale)) pars$scale <- softplus(pars$scale)
  }
  if (dname == "lnorm" && !is.null(pars$sdlog)) {
    pars$sdlog <- softplus(pars$sdlog)
  }
  if (dname == "pois" && !is.null(pars$lambda)) {
    pars$lambda <- softplus(pars$lambda)
  }
  if (dname %in% c("bern", "binom") && !is.null(pars$prob)) {
    pars$prob <- 1 / (1 + exp(-pars$prob))
  }
  if (dname == "binom" && !is.null(pars$size)) {
    pars$size <- softplus(pars$size)
  }
  if (dname == "beta") {
    if (!is.null(pars$shape1)) pars$shape1 <- softplus(pars$shape1)
    if (!is.null(pars$shape2)) pars$shape2 <- softplus(pars$shape2)
  }
  if (dname == "logis" && !is.null(pars$scale)) {
    pars$scale <- softplus(pars$scale)
  }
  pars
}

safe_support <- function(x, dname, pars = list()) {
  valid <- switch(dname,
    exp     = x > 0,
    gamma   = x > 0,
    weibull = x > 0,
    lnorm   = x > 0,
    beta    = (x > 0) & (x < 1),
    pois    = (x >= 0) & (abs(x - round(x)) < 1e-8),
    bern    = x %in% c(0, 1),
    binom   = {
      size <- if (!is.null(pars$size)) pars$size else Inf
      (x >= 0) & (x <= size) & (abs(x - round(x)) < 1e-8)
    },
    rep(TRUE, length(x))
  )
  if (!all(valid)) stop("value outside support")
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

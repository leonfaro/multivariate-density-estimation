
# Triangular transport utilities
#
# - **Input:** `config` list with `K` components.  Each `config[[k]]` defines
#   the conditional distribution `distr` of \(X_k\mid X_{<k}\) and optionally a
#   parameter function `parm`.
# - **Output:** helper functions
#     * `pdf_k`, `cdf_k`, `qtf_k` for density, distribution and quantile of each
#       conditional law.
#     * `S_inv` to compute \(X \sim \pi\) from reference variables \(U_\eta\) or
#       \(Z_\eta\) including log-Jacobian values.
#     * `pi_sample` to draw samples from the target by sequential inversion.
# - **Algorithm:** for \(k=1,\dots,K\)
#     1. Evaluate `config[[k]]$parm` on \(X_{<k}\) to obtain distribution
#        parameters.
#     2. Map \(U_{\eta,k}\) through `qtf_k` to get \(X_k\) and accumulate
#        `logd[k]` via `pdf_k`.
#   The monotonicity of each step ensures the map `S_inv` is invertible.

source("00_setup.R")

safe_pars <- function(pars, ...) pars
safe_support <- function(x, ...) x


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
  xk <- xk
  dens <- do.call(dist_fun("d", cfg[[k]]$distr),
                  c(list(x = xk), pars, list(log = log)))
  dens
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
  res
}

det_J <- function(logd) rowSums(logd)
loglik <- function(Z_eta, logdet) {
  -0.5 * rowSums(Z_eta^2) - (ncol(Z_eta)/2) * log(2 * pi) + logdet
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

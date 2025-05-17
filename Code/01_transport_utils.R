source("00_setup.R")

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
  ck$parm(x_df)
}

pdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  pars <- get_pars(k, x_prev, cfg)
  dens <- do.call(dist_fun("d", cfg[[k]]$distr),
                  c(list(x = xk), pars, list(log = FALSE)))
  if (log) safe_logdens(dens) else dens
}

cdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  dname <- cfg[[k]]$distr
  pars  <- get_pars(k, x_prev, cfg)
  args  <- c(if (dname == "sn") list(x = xk) else list(q = xk),
             pars,
             list(log.p = FALSE))
  cdfv <- do.call(dist_fun("p", dname), args)
  cdfv <- safe_cdf(cdfv)
  if (log) log(cdfv) else cdfv
}

qtf_k <- function(k, u, x_prev, cfg, eps = EPS) {
  u <- pmin(1 - eps, pmax(eps, u))
  dname <- cfg[[k]]$distr
  pars <- get_pars(k, x_prev, cfg)
  if (dname == "sn") {
    res <- tryCatch(
      do.call(dist_fun("q", dname),
              list(p = u, dp = c(unlist(pars), tau = 0), solver = "RFB")),
      error = function(e)
        do.call(dist_fun("q", dname),
                list(p = u, dp = c(unlist(pars), tau = 0), solver = "NR"))
    )
  } else {
    res <- do.call(dist_fun("q", dname), c(list(p = u), pars))
  }
  res[is.infinite(res)] <- sign(res[is.infinite(res)]) * 1e6
  res[is.na(res)] <- 0
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
S_inv <- function(U_eta, cfg = config) {
  U_eta <- apply(U_eta, 2, function(u) pmin(1 - EPS, pmax(EPS, u)))
  n <- nrow(U_eta)
  X_pi <- matrix(NA_real_, n, K)
  logd <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    x_prev <- numeric(0)
    for (k in seq_len(K)) x_prev[k] <- qtf_k(k, U_eta[i, k], x_prev, cfg)
    X_pi[i, ] <- x_prev
    for (k in seq_len(K)) {
      logd[i, k] <- pdf_k(k, X_pi[i, k], X_pi[i, seq_len(k-1)], cfg, log = TRUE) -
                    dnorm(qnorm(U_eta[i, k]), log = TRUE)
    }
  }
  list(X_pi = X_pi, U_eta = U_eta, Z_eta = qnorm(U_eta), logd = logd)
}
pi_sample <- function(N, cfg = config) {
  ref <- eta_sample(N)
  inv <- S_inv(ref$U_eta, cfg)
  list(X_pi = inv$X_pi, U_eta = ref$U_eta,
       Z_eta = inv$Z_eta, logd = inv$logd)
}

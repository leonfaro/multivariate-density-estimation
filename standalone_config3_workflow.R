# Standalone workflow using configuration lists
# Provides functions for sampling, density evaluation, and triangular mapping

config3 <- list(
  list(distr = "norm",  parm = NULL),
  list(distr = "exp",   parm = function(d) list(rate = abs(d$X1) + 1e-8)),
  list(distr = "gamma", parm = function(d) list(shape = d$X2, rate = 1))
)

get_funcs <- function(distr) {
  switch(distr,
    norm  = list(r = rnorm,  d = dnorm,  p = pnorm,  q = qnorm),
    exp   = list(r = rexp,   d = dexp,   p = pexp,   q = qexp),
    gamma = list(r = rgamma, d = dgamma, p = pgamma, q = qgamma),
    stop("unsupported distribution")
  )
}

sample_eta <- function(n, K) {
  matrix(rnorm(n * K), nrow = n, ncol = K,
         dimnames = list(NULL, paste0("Z", seq_len(K))))
}

eta_density <- function(z) {
  apply(dnorm(z), 1, prod)
}

sample_pi <- function(n, config) {
  K <- length(config)
  X <- matrix(NA_real_, nrow = n, ncol = K,
              dimnames = list(NULL, paste0("X", seq_len(K))))
  for (i in seq_len(n)) {
    for (k in seq_len(K)) {
      if (k == 1) {
        d_prev <- data.frame()
      } else {
        d_prev <- as.data.frame(t(X[i, 1:(k - 1)]))
        colnames(d_prev) <- paste0("X", seq_len(k - 1))
      }
      parms <- if (is.null(config[[k]]$parm)) list() else config[[k]]$parm(d_prev)
      rfun <- get_funcs(config[[k]]$distr)$r
      X[i, k] <- do.call(rfun, c(list(1), parms))
    }
  }
  X
}

S <- function(x, config) {
  K <- length(config)
  z <- matrix(NA_real_, nrow = nrow(x), ncol = K,
              dimnames = list(NULL, paste0("Z", seq_len(K))))
  for (k in seq_len(K)) {
    if (k == 1) {
      d_prev <- NULL
    } else {
      d_prev <- as.data.frame(x[, seq_len(k - 1), drop = FALSE])
      colnames(d_prev) <- paste0("X", seq_len(k - 1))
    }
    parms <- if (is.null(config[[k]]$parm)) list() else config[[k]]$parm(d_prev)
    cdf_fun <- get_funcs(config[[k]]$distr)$p
    p_val <- do.call(cdf_fun, c(list(x[, k]), parms))
    z[, k] <- qnorm(p_val)
  }
  z
}

S_inv <- function(z, config) {
  K <- length(config)
  X <- matrix(NA_real_, nrow = nrow(z), ncol = K,
              dimnames = list(NULL, paste0("X", seq_len(K))))
  for (i in seq_len(nrow(z))) {
    for (k in seq_len(K)) {
      if (k == 1) {
        d_prev <- data.frame()
      } else {
        d_prev <- as.data.frame(t(X[i, 1:(k - 1)]))
        colnames(d_prev) <- paste0("X", seq_len(k - 1))
      }
      parms <- if (is.null(config[[k]]$parm)) list() else config[[k]]$parm(d_prev)
      qfun <- get_funcs(config[[k]]$distr)$q
      X[i, k] <- do.call(qfun, c(list(pnorm(z[i, k])), parms))
    }
  }
  X
}

pi_density <- function(x, config) {
  z <- S(x, config)
  phi_prod <- eta_density(z)
  K <- length(config)
  jac <- matrix(NA_real_, nrow = nrow(x), ncol = K)
  for (k in seq_len(K)) {
    if (k == 1) {
      d_prev <- NULL
    } else {
      d_prev <- as.data.frame(x[, seq_len(k - 1), drop = FALSE])
      colnames(d_prev) <- paste0("X", seq_len(k - 1))
    }
    parms <- if (is.null(config[[k]]$parm)) list() else config[[k]]$parm(d_prev)
    d_fun <- get_funcs(config[[k]]$distr)$d
    f_val <- do.call(d_fun, c(list(x[, k]), parms))
    jac[, k] <- f_val / dnorm(z[, k])
  }
  phi_prod * apply(jac, 1, prod)
}

sample_conditional <- function(n, config, x_fixed = rep(NA_real_, length(config))) {
  Z <- sample_eta(n, length(config))
  X <- matrix(NA_real_, nrow = n, ncol = length(config),
              dimnames = list(NULL, paste0("X", seq_len(length(config)))))
  for (i in seq_len(n)) {
    for (k in seq_len(length(config))) {
      if (!is.na(x_fixed[k])) {
        X[i, k] <- x_fixed[k]
        next
      }
      if (k == 1) {
        d_prev <- data.frame()
      } else {
        d_prev <- as.data.frame(t(X[i, 1:(k - 1)]))
        colnames(d_prev) <- paste0("X", seq_len(k - 1))
      }
      parms <- if (is.null(config[[k]]$parm)) list() else config[[k]]$parm(d_prev)
      qfun <- get_funcs(config[[k]]$distr)$q
      X[i, k] <- do.call(qfun, c(list(pnorm(Z[i, k])), parms))
    }
  }
  X
}

# Example usage ------------------------------------------------------------
set.seed(42)
K <- length(config3)
X_train <- sample_pi(50, config3)
X_test  <- sample_pi(50, config3)
stopifnot(!any(is.na(X_train)))
stopifnot(!any(is.na(X_test)))
pi_density(matrix(X_train[1, ], nrow = 1), config3)


suppressPackageStartupMessages({
  library(extraDistr)
  library(sn)
  library(trtf)
})


config25 <- list(
  list(distr = "norm",
       parm  = NULL),

  list(distr = "t",
       parm  = \(d) list(df = 3 + 0.5*d$X1)),

  list(distr = "laplace",
       parm  = \(d) list(m = 0.3*d$X2,
                         s = clip(1 + 0.1*d$X2))),

  list(distr = "logis",
       parm  = \(d) list(location = 0.2*d$X3,
                         scale    = clip(1 + 0.05*d$X3))),

  list(distr = "cauchy",
       parm  = \(d) list(location = 0.1*d$X4,
                         scale    = 0.5 + 0.05*d$X4)),

  list(distr = "exp",
       parm  = \(d) list(rate = clip(1 + 0.1*d$X5))),

  list(distr = "gamma",
       parm  = \(d) list(shape = clip(2 + 0.2*d$X6),
                         rate  = 1)),

  list(distr = "weibull",
       parm  = \(d) list(shape = clip(2 + 0.2*d$X7),
                         scale = clip(1 + 0.1*d$X7))),

  list(distr = "lnorm",
       parm  = \(d) list(meanlog = 0.3*d$X8,
                         sdlog   = 0.5 + 0.05*d$X8)),

  list(distr = "chisq",
       parm  = \(d) list(df = clip(4 + 0.2*d$X9))),

  list(distr = "f",
       parm  = \(d) list(df1 = 5,
                         df2 = 6 + 0.3*d$X10)),

  list(distr = "beta",
       parm  = \(d) { s <- 2 + 0.1*d$X11
                     list(shape1 = s, shape2 = s) }),

  list(distr = "beta",
       parm  = \(d) list(shape1 = 0.5 + 0.05*d$X12,
                         shape2 = 5   + 0.2*d$X12)),

  list(distr = "gumbel",
       parm  = \(d) list(mu    = 0.5*d$X13,
                         sigma = clip(1 + 0.05*d$X13))),

  list(distr = "gpd",
       parm  = \(d) list(mu = 0,
                         sigma = 1 + 0.05*d$X14,
                         xi = 0.1)),

  list(distr = "norm",
       parm  = \(d) list(mean = 0.2*d$X15,
                         sd   = 1)),

  list(distr = "power",
       parm  = \(d) list(alpha = clip(1.5 + 0.1*d$X16),
                         beta  = clip(1))),

  list(distr = "invgamma",
       parm  = \(d) list(alpha = clip(3 + 0.1*d$X17),
                         beta  = clip(1))),

  list(distr = "lomax",
       parm  = \(d) list(lambda = clip(1 + 0.05*d$X18),
                         kappa  = clip(2))),

  list(distr = "rayleigh",
       parm  = \(d) list(sigma = clip(1 + 0.05*d$X19))),

  list(distr = "tnorm",
       parm  = \(d) list(mean = 0.1*d$X20, sd = 1, a = 0, b = Inf)),

  list(distr = "norm",
       parm  = \(d) list(mean = 0.2*d$X21,
                         sd   = 1)),

  list(distr = "hnorm",
       parm  = \(d) list(sigma = clip(1 + 0.05*d$X22))),

  list(distr = "norm",
       parm  = \(d) list(mean = 0.2*d$X23, sd = 1)),

  list(distr = "norm",
       parm  = \(d) list(mean = 0.2*d$X24, sd = 1))
)



config <- config25
EPS <- 1e-8
clip <- function(x, lo = 1e-6, hi = 1e6) pmin(hi, pmax(lo, x))


dist_fun <- function(pref, d) get(paste0(pref, d))


get_pars <- function(k, x_prev, cfg) {
  ck <- cfg[[k]]
  if (is.null(ck$parm)) return(list())
  names(x_prev) <- paste0("X", seq_along(x_prev))
  ck$parm(as.data.frame(as.list(x_prev)))
}


pdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  do.call(dist_fun("d", cfg[[k]]$distr),
          c(list(x = xk), get_pars(k, x_prev, cfg), list(log = log)))
}


cdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  dname <- cfg[[k]]$distr
  fun <- dist_fun("p", dname)
  pars <- get_pars(k, x_prev, cfg)
  if (dname == "sn") {
    args <- c(list(x = xk), pars, list(log.p = log))
  } else {
    args <- c(list(q = xk), pars, list(log.p = log))
  }
  do.call(fun, args)
}


qtf_k <- function(k, u, x_prev, cfg, eps = EPS) {
  u <- pmin(1 - eps, pmax(eps, u))
  dname <- cfg[[k]]$distr
  pars <- get_pars(k, x_prev, cfg)
  if (dname == "sn") {
    tryCatch(
      do.call(dist_fun("q", dname),
              list(p = u,
                   dp = c(unlist(pars), tau = 0),
                   solver = "RFB")),
      error = function(e) {
        do.call(dist_fun("q", dname),
                list(p = u,
                     dp = c(unlist(pars), tau = 0),
                     solver = "NR"))
      }
    )
  } else {
    do.call(dist_fun("q", dname),
            c(list(p = u), pars))
  }
}


pi_density <- function(x, cfg = config, log = FALSE) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  lp <- apply(x, 1, function(row) {
    sum(vapply(seq_along(row),
               function(k)
                 pdf_k(k, row[k], row[seq_len(k-1)], cfg, log = TRUE),
               numeric(1)))
  })
  if (log) lp else exp(lp)
}


eta_sample <- function(N = 1L, K = length(config)) {
  Z <- matrix(rnorm(N*K), nrow = N)
  U <- pnorm(Z)
  list(U = U, Z = Z)
}


S_map <- function(X, cfg = config, eps = EPS) {
  if (is.null(dim(X))) X <- matrix(X, nrow = 1)
  logU <- t(apply(X, 1, function(row) {
    vapply(seq_along(row), function(k)
      cdf_k(k, row[k], row[seq_len(k-1)], cfg, log = TRUE),
      numeric(1))
  }))
  logU[logU < log(eps)] <- log(eps)
  logU[logU > log(1 - eps)] <- log(1 - eps)
  Z <- qnorm(logU, log.p = TRUE)
  U <- exp(logU)
  list(U = U, Z = Z)
}


S_inv <- function(U, cfg = config, eps = EPS) {
  if (is.null(dim(U))) U <- matrix(U, nrow = 1)
  U[U < eps] <- eps
  U[U > 1 - eps] <- 1 - eps
  Z <- qnorm(U); n <- nrow(U); K <- ncol(U)
  X <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    x_prev <- numeric(0)
    for (k in seq_len(K)) {
      x_prev[k] <- qtf_k(k, U[i,k], x_prev[seq_len(k-1)], cfg)
    }
    X[i, ] <- x_prev
  }
  logd <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    for (k in seq_len(K)) {
      logd[i,k] <-
        pdf_k(k, X[i,k], X[i, seq_len(k-1)], cfg, log = TRUE) -
        dnorm(Z[i,k], log = TRUE)
    }
  }
  list(X = X, Z = Z, logd = logd)
}


det_J <- function(logd) rowSums(logd)
loglik <- function(Z, logdet) {
  -0.5 * rowSums(Z^2) - (ncol(Z)/2)*log(2*pi) + logdet
}


round_trip_test <- function(U, cfg = config, tol = 1e-8) {
  inv <- S_inv(U, cfg)
  fwd <- S_map(inv$X, cfg)
  c(maxZ = max(abs(fwd$Z - inv$Z)),
    maxX = max(abs(S_inv(fwd$U, cfg)$X - inv$X))) <= tol
}


log_pi <- function(x, cfg = config) pi_density(x, cfg, log = TRUE)


pi_sample <- function(N = 1L, cfg = config) {
  ref <- eta_sample(N, length(cfg))
  inv <- S_inv(ref$U, cfg)
  list(X = inv$X, U = ref$U, Z = inv$Z, logd = inv$logd)
}


sample_and_loglik <- function(N = 1L, cfg = config) {
  samp <- pi_sample(N, cfg)
  detJ <- det_J(samp$logd)
  ll <- loglik(samp$Z, detJ)
  cbind(samp$X, loglik = ll)
}


cond_sample <- function(x_fixed, k, N = 1L, cfg = config, eps = EPS) {
  K <- length(cfg)
  stopifnot(k >= 1, k < K, length(x_fixed) == k)
  U_fixed <- numeric(k)
  for (j in seq_len(k)) {
    p <- exp(cdf_k(j, x_fixed[j], x_fixed[seq_len(j-1)], cfg, log = TRUE))
    U_fixed[j] <- pmin(1 - eps, pmax(eps, p))
  }
  Z_fixed <- qnorm(U_fixed)
  X_out <- matrix(NA_real_, N, K)
  for (i in seq_len(N)) {
    z <- rnorm(K); z[1:k] <- Z_fixed
    u <- pnorm(z)
    X_out[i,] <- S_inv(matrix(u, nrow=1), cfg)$X
  }
  colnames(X_out) <- paste0("X", 1:K)
  X_out
}


set.seed(2044)
N <- 50
K <- length(config)


samp <- pi_sample(N, config)


rt <- round_trip_test(samp$U, config)
cat('Round-trip test (maxZ,maxX <= tol):', rt, '\n')


cat('Summary of X matrix:\n')
print(summary(as.data.frame(samp$X)))

cat('Summary of U matrix (flattened):\n')
print(summary(as.vector(samp$U)))
cat('Summary of Z matrix (flattened):\n')
print(summary(as.vector(samp$Z)))


logd_vals <- as.vector(samp$logd)
cat('Summary of component log-densities:\n')
print(summary(logd_vals))
detJ_vals <- det_J(samp$logd)
cat('Summary of log-Jacobian det_J:\n')
print(summary(detJ_vals))


n_show <- 10
logd_show <- samp$logd[1:n_show, ]
ll_vals <- loglik(samp$Z, detJ_vals)
cat('Summary of log-likelihoods:\n')
print(summary(ll_vals))


cat('Distribution types for each dimension:\n')
print(data.frame(k = 1:K, distr = sapply(config, `[[`, 'distr')))


X1 <- samp$X[1,]
U1 <- samp$U[1,]
Z1 <- samp$Z[1,]

x_back <- S_inv(matrix(U1, nrow=1), config)$X
cat('Max abs diff for X back-transform (obs1):', max(abs(x_back - X1)), '\n')


cdf_vals <- vapply(seq_len(K), function(k) {
  p <- exp(cdf_k(k, X1[k], X1[seq_len(k-1)], config, log = TRUE))
  pmin(1 - EPS, pmax(EPS, p))
}, numeric(1))
cat('Summary of conditional CDFs for obs1:\n')
print(summary(cdf_vals))


q_back1 <- vapply(seq_len(K), function(k) {
  qtf_k(k, cdf_vals[k], X1[seq_len(k-1)], config)
}, numeric(1))
max_diff <- max(abs(q_back1 - X1))
cat('Max abs diff x from quantile back-transform (obs1):', max_diff, '\n')

diff_vals <- abs(q_back1 - X1)
k_max <- which.max(diff_vals)
cat('Dimension mit max Abweichung:', k_max, 'Diff:', diff_vals[k_max], '\n')


cat('First 3 rows of Z matrix:\n')
print(head(samp$Z, 3))


cat('High-dimensional EDA done.\n')
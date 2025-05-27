library(testthat)

set.seed(123)
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = d$X1)),
  list(distr = "beta", parm = function(d) list(shape1 = d$X2, shape2 = 1))
)
source("../../00_setup.R", chdir = TRUE, local = TRUE)
source("../../02_sampling.R", chdir = TRUE, local = TRUE)

forward_S <- function(X, cfg) {
  n <- nrow(X)
  Z <- matrix(NA_real_, n, K)
  logd <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    x_prev <- numeric(0)
    for (k in seq_len(K)) {
      logu <- cdf_k(k, X[i, k], x_prev, cfg, log = TRUE)
      Z[i, k] <- qnorm(logu, log.p = TRUE)
      logd[i, k] <- pdf_k(k, X[i, k], x_prev, cfg, log = TRUE) -
        dnorm(Z[i, k], log = TRUE)
      x_prev <- c(x_prev, X[i, k])
    }
  }
  list(Z = Z, logd = logd)
}

test_that("push-forward preserves normality", {
  samp <- pi_sample(500, config)
  fwd <- forward_S(samp$X_pi, config)
  for (k in seq_len(K)) {
    ks <- suppressWarnings(ks.test(fwd$Z[, k], "pnorm"))
    expect_gt(ks$p.value, 0.05)
  }
})


test_that("change of variables identity holds", {
  samp <- pi_sample(400, config)
  fwd <- forward_S(samp$X_pi, config)
  logdet <- logdet_J(fwd$logd)
  lhs <- loglik(fwd$Z, logdet)
  true_ll <- rowSums(sapply(seq_len(K), function(k)
    pdf_k(k, samp$X_pi[, k],
          if (k > 1) samp$X_pi[, 1:(k - 1), drop = FALSE] else numeric(0),
          config, log = TRUE)))
  expect_true(abs(mean(lhs - true_ll)) / max(1, abs(mean(true_ll))) < 0.02)
})


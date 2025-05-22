library(testthat)
source("../../00_setup.R", chdir = TRUE)

# Extreme X_k values should not produce NaN log-densities

test_that("pdf_k handles extreme X_k", {
  extreme <- 1e10
  v1 <- pdf_k(1, extreme, numeric(0), config, log = TRUE)
  v2 <- pdf_k(2, extreme, extreme, config, log = TRUE)
  v3 <- pdf_k(3, extreme, c(extreme, extreme), config, log = TRUE)
  expect_false(any(is.na(c(v1, v2, v3))))
})

# Very negative log-Jacobians still yield finite log-likelihood

test_that("loglik finite for very negative log-Jacobians", {
  z <- matrix(0, nrow = 1, ncol = K)
  ld <- matrix(-1000, nrow = 1, ncol = K)
  ll <- loglik(z, logdet_J(ld))
  expect_true(is.finite(ll))
  expect_lt(ll, -900)
})

# NaN entries in the log-Jacobian propagate to the likelihood

test_that("NaN log-Jacobians propagate", {
  z <- matrix(0, nrow = 1, ncol = K)
  ld <- matrix(0, nrow = 1, ncol = K)
  ld[1, 1] <- NaN
  ll <- loglik(z, logdet_J(ld))
  expect_true(is.nan(ll))
})

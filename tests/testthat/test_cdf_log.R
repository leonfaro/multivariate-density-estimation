library(testthat)
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp", parm = function(d) list(rate = d$X1))
)
source("../../00_setup.R", chdir = TRUE, local = TRUE)

# check log output

test_that("cdf_k uses log.p when supported", {
  val <- 0.2
  expect_equal(cdf_k(1, val, numeric(0), config, log = TRUE),
               pnorm(val, log.p = TRUE))
  rate <- softplus(0.5)
  expect_equal(cdf_k(2, val, 0.5, config, log = TRUE),
               pexp(val, rate = rate, log.p = TRUE))
})

test_that("cdf_k returns probabilities when log = FALSE", {
  val <- 0.3
  expect_equal(cdf_k(1, val, numeric(0), config, log = FALSE),
               pnorm(val, log.p = FALSE))
})

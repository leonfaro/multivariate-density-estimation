library(testthat)
source("../../standalone_config3_workflow.R")

set.seed(1)
X <- sample_pi(10, config3)

test_that("sample_pi dimension", {
  expect_equal(dim(X), c(10, 3))
  expect_false(any(is.na(X)))
})

d <- pi_density(matrix(X[1,], nrow=1), config3)

test_that("pi_density positive", {
  expect_true(is.finite(d))
  expect_gt(d, 0)
})

config4 <- c(
  config3,
  list(list(distr = "norm", parm = function(d) list(mean = d$X3, sd = 1)))
)

set.seed(2)
X4 <- sample_pi(5, config4)
d4 <- pi_density(matrix(X4[1,], nrow = 1), config4)

test_that("config4 works", {
  expect_equal(dim(X4), c(5, 4))
  expect_true(is.finite(d4))
  expect_gt(d4, 0)
})

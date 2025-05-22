library(testthat)
source("../../Prototype/standalone_config3_workflow.R")

set.seed(1)
X <- sample_pi(10)

test_that("sample_pi dimension", {
  expect_equal(dim(X), c(10, 3))
  expect_false(any(is.na(X)))
})

d <- pi_density(matrix(X[1,], nrow=1))

test_that("pi_density positive", {
  expect_true(is.finite(d))
  expect_gt(d, 0)
})

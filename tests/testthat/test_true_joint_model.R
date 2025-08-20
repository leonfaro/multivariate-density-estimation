source("helper_config.R")
source("../../models/true_joint_model.R")

set.seed(1)
X <- matrix(runif(40, 0.1, 0.9), ncol = length(config))

test_that("true_joint_logdensity_by_dim returns finite matrix", {
  logmat <- true_joint_logdensity_by_dim(config, X)
  expect_equal(dim(logmat), dim(X))
  expect_true(all(is.finite(logmat)))
})

test_that("logL_TRUE_JOINT_dim matches column means", {
  logmat <- true_joint_logdensity_by_dim(config, X)
  val <- logL_TRUE_JOINT_dim(config, X)
  expect_equal(val, -colMeans(logmat), tolerance = 1e-12)
})

test_that("gamma shape1/shape2 mapping and clamps work", {
  cfg <- list(list(distr = "gamma", parm = function(prev) list(shape1 = -1, shape2 = -Inf)))
  Xbad <- matrix(-5, nrow = 3, ncol = 1)
  ll <- true_joint_logdensity_by_dim(cfg, Xbad)
  expect_true(all(is.finite(ll)))
})


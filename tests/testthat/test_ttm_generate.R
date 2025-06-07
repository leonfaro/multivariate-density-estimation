source("../../00_globals.R")
source("../../01_data_generation.R")

cfg <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp", parm = function(d) list(rate = softplus(d$X1)))
)

set.seed(123)
res1 <- TTM_generate(cfg, N = 30, seed = 123)
set.seed(123)
res2 <- TTM_generate(cfg, N = 30, seed = 123)

z1 <- t(apply(res1$X, 1L, S_forward, theta = res1$theta_hat))

z_rand <- matrix(rnorm(200), ncol = 2)
x_round <- t(apply(z_rand, 1L, R_inverse, theta = res1$theta_hat))
z_back <- t(apply(x_round, 1L, S_forward, theta = res1$theta_hat))

test_that("reproducibility", {
  expect_equal(res1$X, res2$X)
})

test_that("map transforms to approx standard normal", {
  expect_true(max(abs(colMeans(z1))) < 0.1)
  expect_true(max(abs(cov(z1) - diag(2))) < 0.1)
})

test_that("round-trip property", {
  expect_true(max(abs(z_back - z_rand)) < 1e-12)
})

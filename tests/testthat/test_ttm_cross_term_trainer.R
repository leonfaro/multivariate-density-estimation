source("helper_config.R")
source("../../models/ttm_separable.R")
bg_old <- basis_g
source("../../models/ttm_cross_term.R")

set.seed(12)
n <- 50
prep <- prepare_data(n, config)
fit <- trainCrossTermMap(prep$S)

(test_that("no name collision with basis_g", {
  expect_true(exists("basis_g"))
  expect_identical(body(basis_g), body(bg_old))
  expect_true(exists(".basis_g_ct"))
}))

mu_tr <- colMeans(prep$S$X_tr)
sigma_tr <- apply(prep$S$X_tr, 2, sd) + .Machine$double.eps

(test_that("trainCrossTermMap returns valid metrics", {
  expect_s3_class(fit$S, "ttm_cross_term")
  expect_equal(fit$S$mu, mu_tr)
  expect_equal(fit$S$sigma, sigma_tr)
  expect_true(fit$time_train > 0)
  expect_true(is.finite(fit$NLL_train))
  expect_true(is.finite(fit$NLL_val))
  expect_true(is.finite(fit$NLL_test))
  expect_true(is.finite(fit$stderr_test))
}))

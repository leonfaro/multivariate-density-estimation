source("helper_config.R")
source("../../04_evaluation.R")

set.seed(11)
prep <- prepare_data(40, config)
fit <- trainMarginalMap(prep$S)

mu_tr <- colMeans(prep$S$X_tr)
sigma_tr <- apply(prep$S$X_tr, 2, sd) + .Machine$double.eps

test_that("trainMarginalMap nutzt nur Trainingsstandardisierung", {
  expect_s3_class(fit$S, "ttm_marginal")
  expect_equal(fit$S$mu, mu_tr)
  expect_equal(fit$S$sigma, sigma_tr)
  expect_true(fit$time_train >= 0)
  expect_true(is.finite(fit$NLL_train))
  expect_true(is.finite(fit$NLL_val))
  expect_true(is.finite(fit$NLL_test))
  expect_true(is.finite(fit$stderr_test))
  expect_true(all(is.finite(fit$S$coeffA)))
})

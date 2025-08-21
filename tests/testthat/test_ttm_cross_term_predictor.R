source("helper_config.R")
source("../../04_evaluation.R")

set.seed(13)
n <- 50
prep <- prepare_data(n, config)
fit <- trainCrossTermMap(prep$S)

(test_that("predict.ttm_cross_term shapes and sums", {
  LD <- predict(fit$S, prep$S$X_te, "logdensity_by_dim")
  LD_joint <- predict(fit$S, prep$S$X_te, "logdensity")
  expect_equal(dim(LD), dim(prep$S$X_te))
  expect_true(max(abs(rowSums(LD) - LD_joint)) < 1e-12)
  expect_true(all(is.finite(LD)))
}))

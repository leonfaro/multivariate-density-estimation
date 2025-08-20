source("helper_config.R")
source("../../04_evaluation.R")

set.seed(16)
n <- 50
prep <- prepare_data(n, config)

set.seed(101)
fit1 <- trainCrossTermMap(prep$S)
set.seed(101)
fit2 <- trainCrossTermMap(prep$S)

(test_that("trainCrossTermMap reproducible", {
  expect_equal(fit1$NLL_test, fit2$NLL_test)
  p1 <- predict(fit1$S, prep$S$X_te, "logdensity")
  p2 <- predict(fit2$S, prep$S$X_te, "logdensity")
  expect_equal(p1, p2, tolerance = 1e-12)
}))

source("helper_config.R")
source("../../04_evaluation.R")

set.seed(15)
prep <- prepare_data(30, config)
fit <- trainCrossTermMap(prep$S)

(test_that("forward KL identity holds", {
  LD <- predict(fit$S, prep$S$X_te, "logdensity_by_dim")
  lhs <- .forwardKLLoss_ct(fit$S, prep$S$X_te)
  rhs <- mean(-rowSums(LD))
  expect_lt(abs(lhs - rhs), 1e-10)
}))

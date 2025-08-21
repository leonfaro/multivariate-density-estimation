source("helper_config.R")
source("../../04_evaluation.R")

set.seed(15)
n <- 50
prep <- prepare_data(n, config)
fit <- trainCrossTermMap(prep$S)

(test_that("forward KL identity holds", {
  LD <- predict(fit$S, prep$S$X_te, "logdensity_by_dim")
  lhs <- .forwardKLLoss_ct(fit$S, prep$S$X_te)
  rhs <- mean(-rowSums(LD) - 0.5 * ncol(prep$S$X_te) * log(2 * pi))
  expect_lt(abs(lhs - rhs), 1e-10)
}))

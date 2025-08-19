source("helper_config.R")
source("../../04_evaluation.R")

set.seed(14)
prep <- prepare_data(30, config)
fit <- trainCrossTermMap(prep$S)
S0 <- fit$S
for (k in seq_along(S0$coeffs)) {
  S0$coeffs[[k]]$alpha[] <- 0
  S0$coeffs[[k]]$beta[] <- 0
}

(test_that("predictor yields constants for zero coeffs", {
  LD0 <- predict(S0, prep$S$X_te, "logdensity_by_dim")
  const <- -0.5 * log(2 * pi) - log(S0$sigma)
  expected <- matrix(rep(const, each = nrow(LD0)), nrow(LD0), byrow = FALSE)
  expect_equal(LD0, expected, tolerance = 1e-12)
}))

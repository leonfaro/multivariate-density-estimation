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
  Xs <- .standardize(S0, prep$S$X_te)
  const <- -0.5 * log(2 * pi)
  expected <- (-0.5) * (Xs^2) + const +
    matrix(-log(S0$sigma), nrow(Xs), ncol(Xs), byrow = TRUE)
  expect_equal(LD0, expected, tolerance = 1e-12, check.attributes = FALSE)
}))

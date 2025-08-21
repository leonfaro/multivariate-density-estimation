source("helper_config.R")
source("../../models/ttm_separable.R")
source("../../models/ttm_cross_term.R")

set.seed(2024)
prep <- prepare_data(100, config)

fit_zero <- trainCrossTermMap(prep$S)
fit_ws   <- trainCrossTermMap(prep$S, warmstart_from_separable = TRUE)

(test_that("warm start converges and improves NLL", {
  conv <- vapply(fit_ws$S$coeffs, function(c) c$convergence, numeric(1))
  expect_true(all(conv == 0))
  expect_true(fit_ws$NLL_val <= fit_zero$NLL_val + 1e-4)
}))

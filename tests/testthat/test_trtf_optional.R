context("TRTF Predict-API invariants (optional, if packages installed)")

test_that("TRTF by-dim and joint agree and are finite", {
  if (!requireNamespace("trtf", quietly = TRUE) || !requireNamespace("tram", quietly = TRUE)) {
    skip("trtf/tram not available in test environment")
  }
  set.seed(101)
  n <- 80L
  prep <- prepare_data(n, config, seed = 101)
  S <- prep$S
  mod <- fit_TRTF(S, config, seed = 101)
  LD <- predict(mod, S$X_te, type = "logdensity_by_dim", cores = 1L, trace = FALSE)
  expect_true(is.matrix(LD))
  expect_equal(dim(LD), dim(S$X_te))
  expect_true(all(is.finite(LD)))
  joint <- predict(mod, S$X_te, type = "logdensity", cores = 1L, trace = FALSE)
  expect_equal(length(joint), nrow(S$X_te))
  expect_lte(max(abs(rowSums(LD) - joint)), 1e-10)
})


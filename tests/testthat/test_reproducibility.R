context("Reproducibility with fixed seeds")

test_that("Marginal and Separable fits reproduce NLL under same seed", {
  set.seed(123)
  n <- 80L
  prep <- prepare_data(n, config, seed = 123)
  S <- prep$S

  # Marginal
  m1 <- trainMarginalMap(S, seed = 7)
  m2 <- trainMarginalMap(S, seed = 7)
  expect_equal(m1$NLL_test, m2$NLL_test, tolerance = 1e-12)

  # Separable (deterministic convex objective)
  s1 <- trainSeparableMap(S, degree_g = 2, lambda = 1e-3, seed = 11)
  s2 <- trainSeparableMap(S, degree_g = 2, lambda = 1e-3, seed = 11)
  expect_equal(s1$NLL_test, s2$NLL_test, tolerance = 1e-10)
})


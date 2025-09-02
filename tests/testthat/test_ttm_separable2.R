context("Separable TTM via core: invariants and deg=0 equivalence to marginal")

test_that("ldb_dim invariants; deg=0 matches marginal within 1e-10", {
  source(file.path(root_path, "models/ttm/ttm_marginal.R"))
  source(file.path(root_path, "models/ttm/ttm_separable.R"))
  set.seed(77)
  N <- 1000L; K <- 3L
  X <- matrix(rnorm(N * K), ncol = K)
  S <- split_data(X, seed = 77)
  fit_m <- fit_ttm(S, algo = "marginal", seed = 77)
  fit_s0 <- fit_ttm(S, algo = "separable", degree_g = 0L, lambda = 0.0, seed = 77)

  LDm <- predict_ttm(fit_m,  S$X_te, type = "logdensity_by_dim")
  LDs <- predict_ttm(fit_s0, S$X_te, type = "logdensity_by_dim")
  expect_true(is.matrix(LDm) && is.matrix(LDs))
  expect_equal(dim(LDm), dim(LDs))
  expect_true(all(is.finite(LDm)))
  expect_true(all(is.finite(LDs)))
  # rowSums property
  expect_lte(max(abs(rowSums(LDs) - predict_ttm(fit_s0, S$X_te, type = "logdensity"))), 1e-10)
  # Equivalence to marginal under deg=0
  expect_lte(max(abs(LDm - LDs)), 1e-10)
})

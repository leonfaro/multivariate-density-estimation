context("New marginal TTM (R/ttm_marginal.R) API, invariants, constants, and Gaussian sanity")

test_that("fit vs fit$S accepted; shapes; idempotency; constants; Gaussian sanity", {
  source(file.path(root_path, "R/ttm_marginal.R"))
  set.seed(123)
  N <- 30000L; K <- 2L
  X <- matrix(rnorm(N * K), ncol = K)
  S <- split_data(X, seed = 123)
  fit <- fit_ttm(S, seed = 123)
  M <- fit$S

  # Both call forms
  LD_by1 <- predict_ttm(fit, S$X_te, type = "logdensity_by_dim")
  LD_by2 <- predict_ttm(M,   S$X_te, type = "logdensity_by_dim")
  expect_equal(LD_by1, LD_by2, tolerance = 1e-12)

  # Additional types and shapes
  Z  <- predict_ttm(M, S$X_te, type = "transform")
  D  <- predict_ttm(M, S$X_te, type = "jac_diag")
  LD <- predict_ttm(M, S$X_te, type = "logdensity_by_dim")
  Lj <- predict_ttm(M, S$X_te, type = "logdensity")
  expect_equal(dim(Z),  dim(S$X_te))
  expect_equal(dim(D),  dim(S$X_te))
  expect_equal(dim(LD), dim(S$X_te))
  expect_equal(length(Lj), nrow(S$X_te))
  expect_true(all(is.finite(Z)))
  expect_true(all(is.finite(D)))
  expect_true(all(is.finite(LD)))
  expect_true(all(is.finite(Lj)))
  expect_lte(max(abs(rowSums(LD) - Lj)), 1e-10)

  # Constants exactly once: rebuild LD from parts
  C <- -0.5 * log(2 * pi)
  LJ <- log(D) - matrix(log(M$sigma), nrow = nrow(D), ncol = ncol(D), byrow = TRUE)
  LD_chk <- (-0.5) * (Z^2) + C + LJ
  expect_equal(LD, LD_chk, tolerance = 1e-12)

  # Idempotency
  LD2 <- predict_ttm(M, S$X_te, type = "logdensity_by_dim")
  expect_equal(LD, LD2, tolerance = 1e-12)

  # Gaussian sanity
  per_dim_nll <- -colMeans(LD)
  target <- 0.5 * log(2 * pi * exp(1))
  expect_true(all(abs(per_dim_nll - target) < 1e-2))
})

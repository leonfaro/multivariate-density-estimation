context("Predict-API invariants for TTM models")

test_that("TTM marginal/separable/cross-term return finite LD with correct shapes and invariants", {
  set.seed(42)
  n <- 60L
  prep <- prepare_data(n, config, seed = 42)
  S <- prep$S

  # Marginal
  fit_m <- trainMarginalMap(S, seed = 42)
  M <- fit_m$S
  LDm <- predict(M, S$X_te, "logdensity_by_dim")
  expect_true(is.matrix(LDm))
  expect_equal(dim(LDm), dim(S$X_te))
  expect_true(all(is.finite(LDm)))
  joint_m <- predict(M, S$X_te, "logdensity")
  expect_equal(length(joint_m), nrow(S$X_te))
  expect_true(all(is.finite(joint_m)))
  expect_lte(max(abs(rowSums(LDm) - joint_m)), 1e-10)

  # Separable
  fit_s <- trainSeparableMap(S, degree_g = 2, lambda = 1e-3, seed = 42)
  Ms <- fit_s$S
  LDs <- predict(Ms, S$X_te, "logdensity_by_dim")
  expect_equal(dim(LDs), dim(S$X_te))
  expect_true(all(is.finite(LDs)))
  joint_s <- predict(Ms, S$X_te, "logdensity")
  expect_true(all(is.finite(joint_s)))
  expect_lte(max(abs(rowSums(LDs) - joint_s)), 1e-10)

  # Cross-term (keep light via small Q)
  oldQ <- getOption("mde.ctm.Q", NULL)
  on.exit({ if (!is.null(oldQ)) options(mde.ctm.Q = oldQ) else options(mde.ctm.Q = NULL) }, add = TRUE)
  options(mde.ctm.Q = 16L)
  fit_c <- trainCrossTermMap(S, degree_g = 2, seed = 42, warmstart_from_separable = TRUE,
                             sep_degree_g = 2, sep_lambda = 1e-3)
  Mc <- fit_c$S
  LDc <- predict(Mc, S$X_te, "logdensity_by_dim")
  expect_equal(dim(LDc), dim(S$X_te))
  expect_true(all(is.finite(LDc)))
  joint_c <- predict(Mc, S$X_te, "logdensity")
  expect_true(all(is.finite(joint_c)))
  expect_lte(max(abs(rowSums(LDc) - joint_c)), 1e-10)
})


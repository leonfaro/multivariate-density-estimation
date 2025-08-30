context("Unified TTM API invariants across {marginal, separable, crossterm}")

test_that("All algos expose types, shapes, constants, rowSums invariants", {
  source(file.path(root_path, "R/ttm_bases.R"))
  source(file.path(root_path, "R/ttm_core.R"))
  source(file.path(root_path, "R/ttm_marginal.R"))
  source(file.path(root_path, "R/ttm_separable.R"))
  source(file.path(root_path, "R/ttm_crossterm.R"))
  set.seed(1)
  N <- 1200L; K <- 3L
  X <- matrix(rnorm(N * K), ncol = K)
  S <- split_data(X, seed = 1)
  algos <- c("marginal", "separable", "crossterm")
  for (algo in algos) {
    if (algo == "marginal") {
      fit <- fit_ttm(S, algo = "marginal", seed = 1)
    } else if (algo == "separable") {
      fit <- fit_ttm(S, algo = "separable", degree_g = 1L, lambda = 1e-3, seed = 1)
    } else if (algo == "crossterm") {
      fit <- fit_ttm(S, algo = "crossterm", deg_g = 0L, df_t = 4L, Q = 8L, lambda = 1e-3, seed = 1)
    }
    Z <- predict_ttm(fit, S$X_te, type = "transform")
    J <- predict_ttm(fit, S$X_te, type = "jac_diag")
    L <- predict_ttm(fit, S$X_te, type = "logdensity_by_dim")
    Lj <- predict_ttm(fit, S$X_te, type = "logdensity")
    expect_equal(dim(Z), dim(S$X_te))
    expect_equal(dim(J), dim(S$X_te))
    expect_equal(dim(L), dim(S$X_te))
    expect_equal(length(Lj), nrow(S$X_te))
    expect_true(all(is.finite(Z)))
    expect_true(all(is.finite(J)))
    expect_true(all(is.finite(L)))
    expect_true(all(is.finite(Lj)))
    expect_lte(max(abs(rowSums(L) - Lj)), 1e-12)
    # Constants exactly once per dimension
    C <- -0.5 * log(2 * pi)
    Sigma <- matrix(log(fit$S$sigma), nrow = nrow(L), ncol = ncol(L), byrow = TRUE)
    L_chk <- (-0.5) * (Z^2) + C + log(J) - Sigma
    expect_equal(L, L_chk, tolerance = 1e-12)
  }
})

test_that("Gaussian marginal NLL per dimension ≈ 1.41894 (±1e-2)", {
  source(file.path(root_path, "R/ttm_bases.R"))
  source(file.path(root_path, "R/ttm_core.R"))
  source(file.path(root_path, "R/ttm_marginal.R"))
  set.seed(2)
  N <- 30000L; K <- 2L
  X <- matrix(rnorm(N * K), ncol = K)
  S <- split_data(X, seed = 2)
  fit <- fit_ttm(S, algo = "marginal", seed = 2)
  L <- predict_ttm(fit, S$X_te, type = "logdensity_by_dim")
  per_dim_nll <- -colMeans(L)
  target <- 0.5 * log(2 * pi * exp(1))
  expect_true(all(abs(per_dim_nll - target) < 1e-2))
})


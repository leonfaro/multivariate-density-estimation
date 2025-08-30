context("Modularity limits: separable→marginal (exact) and crossterm→separable (near)")

test_that("Separable with deg=0 matches Marginal exactly (joint LD)", {
  source(file.path(root_path, "R/ttm_bases.R"))
  source(file.path(root_path, "R/ttm_core.R"))
  source(file.path(root_path, "R/ttm_marginal.R"))
  source(file.path(root_path, "R/ttm_separable.R"))
  set.seed(1)
  N <- 1200L; K <- 3L
  X <- matrix(rnorm(N * K), ncol = K)
  S <- split_data(X, seed = 1)
  fit_m <- fit_ttm(S, algo = "marginal",  seed = 1)
  fit_s <- fit_ttm(S, algo = "separable", degree_g = 0L, lambda = 0.0, seed = 1)
  Lm <- predict_ttm(fit_m, S$X_te, type = "logdensity")
  Ls <- predict_ttm(fit_s, S$X_te, type = "logdensity")
  expect_true(all(is.finite(Lm)))
  expect_true(all(is.finite(Ls)))
  expect_equal(Lm, Ls, tolerance = 1e-10)
})

test_that("Cross-term with h(t) only is close to separable deg=0 (joint LD)", {
  source(file.path(root_path, "R/ttm_bases.R"))
  source(file.path(root_path, "R/ttm_core.R"))
  source(file.path(root_path, "R/ttm_marginal.R"))
  source(file.path(root_path, "R/ttm_separable.R"))
  source(file.path(root_path, "R/ttm_crossterm.R"))
  set.seed(1)
  N <- 600L; K <- 2L
  X <- matrix(rnorm(N * K), ncol = K)
  S <- split_data(X, seed = 1)
  fit_s <- fit_ttm(S, algo = "separable", degree_g = 0L, lambda = 0.0, seed = 1)
  # h independent of x_<k> via deg_g=0; moderate quadrature
  fit_c <- fit_ttm(S, algo = "crossterm", deg_g = 0L, df_t = 6L, Q = 16L, lambda = 1e-3, seed = 1)
  Ls <- predict_ttm(fit_s, S$X_te, type = "logdensity")
  Lc <- predict_ttm(fit_c, S$X_te, type = "logdensity")
  expect_true(all(is.finite(Ls)))
  expect_true(all(is.finite(Lc)))
  # RowSums invariant sanity on cross-term
  Lc_by <- predict_ttm(fit_c, S$X_te, type = "logdensity_by_dim")
  expect_lte(max(abs(rowSums(Lc_by) - Lc)), 1e-12)
  # Near-identical on test
  diff <- mean(abs(Lc - Ls))
  expect_lt(diff, 1e-3)
})


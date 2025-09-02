context("Unified bases: build_f/d_build_f, build_g, build_h, GL nodes")

test_that("build_f and d_build_f shapes and derivative numeric check", {
  source(file.path(root_path, "models/ttm/ttm_bases.R"))
  set.seed(7)
  x <- rnorm(100)
  F <- build_f(x)
  dF <- d_build_f(x)
  expect_true(is.matrix(F) && is.matrix(dF))
  expect_equal(dim(F), c(length(x), 2L))
  expect_equal(dim(dF), c(length(x), 2L))
  # Numeric derivative check
  eps <- 1e-6
  Fp <- build_f(x + eps)
  Fm <- build_f(x - eps)
  num <- (Fp - Fm) / (2 * eps)
  expect_lte(max(abs(num - dF)), 1e-8)
})

test_that("build_g constructs dense polynomial features with correct dims", {
  source(file.path(root_path, "models/ttm/ttm_bases.R"))
  set.seed(3)
  N <- 20L; p <- 3L; deg <- 3L
  Xp <- matrix(rnorm(N * p), ncol = p)
  G <- build_g(Xp, deg = deg)
  expect_true(is.matrix(G) && all(is.finite(G)))
  expect_equal(dim(G), c(N, 1L + p * deg))
})

test_that("build_h returns tensor of B-splines(t) and poly(X_prev)", {
  source(file.path(root_path, "models/ttm/ttm_bases.R"))
  set.seed(9)
  N <- 25L; p <- 2L
  t <- runif(N)
  Xp <- matrix(rnorm(N * p), ncol = p)
  spec <- list(df = 6L, degree = 3L, deg_g = 2L)
  H <- build_h(t, Xp, spec)
  expect_true(is.matrix(H) && all(is.finite(H)))
  # dims: N x (df * (1 + p*deg_g))
  expect_equal(dim(H), c(N, spec$df * (1L + p * spec$deg_g)))
})

test_that("Gauss–Legendre nodes on [0,1] sum weights to 1 and are in (0,1)", {
  source(file.path(root_path, "models/ttm/ttm_bases.R"))
  gl1 <- gauss_legendre_nodes(1L)
  expect_equal(gl1$nodes, 0.5)
  expect_equal(gl1$weights, 1.0)
  gl <- gauss_legendre_nodes(8L)
  expect_true(all(gl$nodes > 0 & gl$nodes < 1))
  expect_lt(abs(sum(gl$weights) - 1), 1e-12)
})

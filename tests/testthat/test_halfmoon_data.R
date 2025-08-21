source("../../scripts/halfmoon_data.R")

context("two-moons dataset")

test_that("generate_two_moons reproducibility and shape", {
  n <- 40L
  noise <- 0.05
  seed <- 7
  X1 <- generate_two_moons(n, noise, seed)
  X2 <- generate_two_moons(n, noise, seed)
  expect_true(is.matrix(X1))
  expect_identical(dim(X1), c(n, 2L))
  expect_identical(colnames(X1), c("x1", "x2"))
  expect_identical(X1, X2)
  expect_false(any(!is.finite(X1)))
})

test_that("make_halfmoon_splits produces consistent splits", {
  n_train <- 60L
  n_test <- 30L
  noise <- 0.1
  seed <- 13
  S1 <- make_halfmoon_splits(n_train, n_test, noise, seed)
  S2 <- make_halfmoon_splits(n_train, n_test, noise, seed)
  expect_identical(S1, S2)
  n_val <- as.integer(max(10, round(0.2 * n_train)))
  expect_identical(dim(S1$X_tr), c(as.integer(n_train - n_val), 2L))
  expect_identical(dim(S1$X_val), c(n_val, 2L))
  expect_identical(dim(S1$X_te), c(n_test, 2L))
  expect_identical(colnames(S1$X_tr), c("x1", "x2"))
  expect_identical(colnames(S1$X_val), c("x1", "x2"))
  expect_identical(colnames(S1$X_te), c("x1", "x2"))
  expect_equal(S1$K, 2L)
  expect_false(any(!is.finite(S1$X_tr)))
  expect_false(any(!is.finite(S1$X_val)))
  expect_false(any(!is.finite(S1$X_te)))
  expect_identical(check_halfmoon_splits(S1), S1)
})

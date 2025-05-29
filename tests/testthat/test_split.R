source("../../02_split.R")

set.seed(1)
X_dummy <- matrix(rnorm(100), nrow = 20)


test_that("train_test_split returns proper dimensions", {
  out <- train_test_split(X_dummy, 0.5, 123)
  expect_equal(nrow(out$X_tr), floor(0.5 * nrow(X_dummy)))
  expect_equal(nrow(out$X_te), nrow(X_dummy) - floor(0.5 * nrow(X_dummy)))
  expect_equal(ncol(out$X_tr), ncol(X_dummy))
  expect_equal(ncol(out$X_te), ncol(X_dummy))
})

test_that("train_test_split uses seed correctly", {
  out1 <- train_test_split(X_dummy, 0.4, 42)
  out2 <- train_test_split(X_dummy, 0.4, 42)
  expect_equal(out1$X_tr, out2$X_tr)
  expect_equal(out1$X_te, out2$X_te)
})

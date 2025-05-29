source("../../00_globals.R")
source("../../01_data_generation.R")
source("../../models/ttm_model.R")

set.seed(42)
G <- setup_global()
X <- gen_samples(G)
n_tr <- floor(G$split_ratio * G$N)
X_tr <- X[seq_len(n_tr), ]
X_te <- X[(n_tr + 1):nrow(X), ]

# basic fit_TTM functionality

test_that("fit_TTM returns expected structure", {
  res <- fit_TTM(X_tr, X_te, G$H_grid[1])
  expect_type(res, "list")
  expect_true(is.numeric(res$theta))
  expect_true(length(res$theta) > 0)
  expect_true(is.numeric(res$h) && length(res$h) == 1)
  expect_true(is.numeric(res$logL_te))
  expect_true(is.finite(res$logL_te))
})

test_that("logL_TTM computes finite value", {
  res <- fit_TTM(X_tr, X_te, G$H_grid[1])
  val <- logL_TTM(res, X_te)
  expect_true(is.finite(val))
})

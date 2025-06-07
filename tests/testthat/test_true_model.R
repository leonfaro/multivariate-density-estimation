source("../../00_globals.R")
source("../../01_data_generation.R")
source("../../models/true_model.R")

set.seed(1)
G <- setup_global()
X <- gen_samples(G)
N_tr <- floor(G$split_ratio * G$n)
X_tr <- X[seq_len(N_tr), ]
X_te <- X[(N_tr + 1):nrow(X), ]

test_that("fit_TRUE returns parameter list", {
  res <- fit_TRUE(X_tr, X_te, G$config)
  expect_type(res, "list")
  expect_true(is.list(res$theta))
  expect_equal(length(res$theta), ncol(X))
  expect_true(is.numeric(res$logL_te))
  expect_true(is.finite(res$logL_te))
})

test_that("logL_TRUE computes finite value", {
  res <- fit_TRUE(X_tr, X_te, G$config)
  val <- logL_TRUE(res, X_te)
  expect_true(is.finite(val))
})

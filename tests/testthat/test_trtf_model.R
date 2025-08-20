source("../../00_globals.R")
source("../../01_data_generation.R")
source("../../02_split.R")
source("../../models/trtf_model.R")

set.seed(1)
n <- 50
X <- Generate_iid_from_config(n, config)
S <- split_data(X, 42)


test_that("fit_TRTF returns valid object", {
  mod <- fit_TRTF(S, config)
  expect_type(mod, "list")
  expect_s3_class(mod, "mytrtf")
  expect_equal(length(mod$ymod), ncol(S$X_tr))
  expect_equal(length(mod$forests), ncol(S$X_tr) - 1)
  expect_true(is.numeric(mod$logL_te))
  expect_true(is.finite(mod$logL_te))
})


test_that("logL_TRTF computes finite value", {
  mod <- fit_TRTF(S, config)
  val <- logL_TRTF(mod, S$X_te)
  expect_true(is.finite(val))
})

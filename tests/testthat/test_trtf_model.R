source("../../00_globals.R")
source("../../01_data_generation.R")
source("../../02_split.R")
source("../../models/trtf_model.R")

set.seed(1)
X <- Generate_iid_from_config(60, config)
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


test_that("predict.mytrtf logdensity aggregation works", {
  mod <- fit_TRTF(S, config, seed = 42)
  ll_dim <- predict(mod, S$X_te, type = "logdensity_by_dim")
  expect_equal(dim(ll_dim), dim(S$X_te))
  expect_true(all(is.finite(ll_dim)))
  ll <- predict(mod, S$X_te, type = "logdensity")
  expect_equal(ll, rowSums(ll_dim), tolerance = 1e-12)
})


test_that("predict.mytrtf ist deterministisch bei seed 42", {
  set.seed(42)
  mod1 <- fit_TRTF(S, config, seed = 42)
  pred1 <- predict(mod1, S$X_te, type = "logdensity")
  set.seed(42)
  mod2 <- fit_TRTF(S, config, seed = 42)
  pred2 <- predict(mod2, S$X_te, type = "logdensity")
  expect_identical(pred1, pred2)
})

source("../../00_globals.R")
source("../../01_data_generation.R")
source("../../02_split.R")
source("../../models/ttm_model.R")

set.seed(1)
G <- setup_global()
G$n <- 40
X <- gen_samples(G)
S <- train_test_split(X, G$split_ratio, G$seed)


test_that("fit_TTM liefert valides Objekt", {
  mod <- fit_TTM(S$X_tr, S$X_te)
  expect_type(mod, "list")
  expect_true(is.numeric(mod$test_logL))
  expect_true(is.finite(mod$test_logL))
})


test_that("logL_TTM gibt endlichen Wert", {
  mod <- fit_TTM(S$X_tr, S$X_te)
  val <- logL_TTM(mod, S$X_te)
  expect_true(is.finite(val))
})

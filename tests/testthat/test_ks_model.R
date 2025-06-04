source("../../00_globals.R")
source("../../01_data_generation.R")
source("../../02_split.R")
source("../../models/ks_model.R")

set.seed(2)
G <- setup_global()
G$N <- 60
X <- gen_samples(G)
S <- train_test_split(X, G$split_ratio, G$seed)


test_that("fit_KS liefert valides Objekt", {
  mod <- fit_KS(S$X_tr, S$X_te, G$config)
  expect_type(mod, "list")
  expect_s3_class(mod, "ks_model")
  expect_true(is.numeric(mod$logL_te))
  expect_true(is.finite(mod$logL_te))
})

test_that("logL_KS gibt endlichen Wert", {
  mod <- fit_KS(S$X_tr, S$X_te, G$config)
  val <- logL_KS(mod, S$X_te)
  expect_true(is.finite(val))
})

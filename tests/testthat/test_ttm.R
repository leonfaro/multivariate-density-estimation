source("../../00_globals.R")
source("../../01_data_generation.R")
source("../../02_split.R")
source("../../models/triangular_transport_map.R")

set.seed(1)
G <- setup_global()
G$N <- 60
X <- gen_samples(G)
S <- train_test_split(X, G$split_ratio, G$seed)


test_that("fit_TTM liefert valides Modell", {
  mod <- fit_TTM(S$X_tr, S$X_te)
  expect_type(mod, "list")
  expect_s3_class(mod, "ttm")
  expect_true(is.numeric(mod$logL_te))
  expect_true(is.finite(mod$logL_te))
  expect_equal(length(mod$mu), ncol(S$X_tr))
})


test_that("logL_TTM berechnet endlichen Wert", {
  mod <- fit_TTM(S$X_tr, S$X_te)
  val <- logL_TTM(mod, S$X_te)
  expect_true(is.finite(val))
})

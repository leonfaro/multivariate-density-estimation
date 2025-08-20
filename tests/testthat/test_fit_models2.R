source("helper_config.R")
source("../../04_evaluation.R")
source("../../models/trtf_model.R")
source("../../models/true_model.R")

set.seed(2)
prep <- prepare_data(30, config)
  mods <- fit_models(prep$S, config)


test_that("fit_models returns list of models", {
  expect_true(is.list(mods$models))
  expect_s3_class(mods$models$trtf, "mytrtf")
   expect_s3_class(mods$models$ttm_cross$S, "ttm_cross_term")
  expect_true(is.list(mods$ll))
  expect_length(mods$ll$true, length(config))
  expect_equal(ncol(mods$ll$ttm_cross), length(config))
  expect_true(all(is.finite(unlist(mods$ll))))
  expect_true(all(mods$times >= 0))
})

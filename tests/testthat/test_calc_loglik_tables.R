source("helper_config.R")
source("../../04_evaluation.R")
source("../../models/trtf_model.R")
source("../../models/ks_model.R")
source("../../models/true_model.R")

set.seed(3)
prep <- prepare_data(30, config)
mods <- fit_models(prep$S, config)

tab <- calc_loglik_tables(mods, config)

test_that("calc_loglik_tables works", {
  expect_true(is.data.frame(tab))
  expect_equal(nrow(tab), length(config) + 1)
  expect_true(all(is.finite(tab$logL_baseline)))
  expect_true(all(is.finite(tab$logL_trtf)))
  expect_true(all(is.finite(tab$logL_ks)))
})

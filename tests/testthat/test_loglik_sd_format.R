source("helper_config.R")
source("../../04_evaluation.R")
source("../../models/trtf_model.R")
source("../../models/ks_model.R")
source("../../models/true_model.R")
source("../../02_split.R")

set.seed(7)
prep <- prepare_data(30, config, c(3,2,1,4))
mods <- fit_models(prep$S, config)

tab_mean <- calc_loglik_tables(mods, config)
tab_sd   <- calc_loglik_sds(mods, prep$S, config)
tab_fmt  <- format_loglik_table(tab_mean, tab_sd)

test_that("calc_loglik_sds works", {
  expect_true(is.data.frame(tab_sd))
  expect_equal(nrow(tab_sd), length(config) + 1)
  expect_true(all(is.finite(tab_sd$logL_trtf)))
})

test_that("format_loglik_table adds strings", {
  expect_true(is.character(tab_fmt$logL_baseline))
  expect_match(tab_fmt$logL_trtf[1], "±")
})

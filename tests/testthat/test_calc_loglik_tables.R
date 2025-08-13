source("helper_config.R")
source("../../04_evaluation.R")
source("../../models/trtf_model.R")
source("../../models/ks_model.R")
source("../../models/true_model.R")
source("../../models/ttm_separable.R")

set.seed(3)
prep <- prepare_data(30, config)
mods <- list(
  true = fit_TRUE(prep$S, config),
  trtf = fit_TRTF(prep$S, config),
  ks   = fit_KS(prep$S, config),
  ttm  = trainMarginalMap(prep$S),
  ttm_sep = trainSeparableMap(prep$S)
)
tab <- calc_loglik_tables(mods, config, prep$S$X_te)

test_that("calc_loglik_tables works", {
  expect_true(is.data.frame(tab))
  expect_equal(nrow(tab), length(config) + 1)
  expect_true(all(grepl("±", tab$`Marginal Map`)))
  expect_true(all(grepl("±", tab$`Separable Map`)))
  expect_false(any(grepl("NA", tab$`Separable Map`)))
})

source("helper_config.R")
source("../../04_evaluation.R")
source("../../models/trtf_model.R")
source("../../models/true_model.R")
source("../../models/ttm_separable.R")

set.seed(7)
n <- 50
prep <- prepare_data(n, config)
  mods <- list(
    true = fit_TRUE(prep$S, config),
    trtf = fit_TRTF(prep$S, config),
    ttm  = trainMarginalMap(prep$S),
    ttm_sep = trainSeparableMap(prep$S)
  )
tab <- calc_loglik_tables(mods, config, prep$S$X_te)

test_that("calc_loglik_tables returns strings with plusminus", {
  expect_true(is.data.frame(tab))
  expect_true(all(grepl("±", tab$`Marginal Map`)))
  expect_true(all(grepl("±", tab$`Separable Map`)))
})

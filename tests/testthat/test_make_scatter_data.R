source("helper_config.R")
source("../../EDA.R")
source("../../models/trtf_model.R")
source("../../models/ks_model.R")
source("../../models/true_model.R")
source("../../02_split.R")

set.seed(4)
prep <- prepare_data(30, config, c(3,2,1,4))
mods <- fit_models(prep$S, config)

scat <- make_scatter_data(mods, prep$S)

test_that("make_scatter_data returns numeric vectors", {
  expect_true(is.list(scat))
  expect_length(scat$ld_base, nrow(prep$S$X_te))
  expect_true(all(is.finite(scat$ld_trtf)))
  expect_true(all(is.finite(scat$ld_ks)))
})

source("helper_config.R")
source("../../04_evaluation.R")

set.seed(1)
res <- prepare_data(20, config, c(3,2,1,4))


test_that("prepare_data returns valid structure", {
  expect_true(is.list(res))
  expect_true(all(c("G", "S", "S_perm", "param_list") %in% names(res)))
  expect_equal(nrow(res$S$X_tr) + nrow(res$S$X_te), 20)
  expect_equal(ncol(res$S$X_tr), length(config))
  expect_equal(ncol(res$S_perm$X_tr), length(config))
})

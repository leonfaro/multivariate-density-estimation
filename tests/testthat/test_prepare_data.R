source("helper_config.R")
source("../../04_evaluation.R")

set.seed(1)
n <- 50
res <- prepare_data(n, config)


test_that("prepare_data returns valid structure", {
  expect_true(is.list(res))
  expect_true(all(c("X", "S") %in% names(res)))
  expect_equal(ncol(res$S$X_tr), length(config))
  expect_equal(nrow(res$S$X_tr) + nrow(res$S$X_val) + nrow(res$S$X_te), n)
})

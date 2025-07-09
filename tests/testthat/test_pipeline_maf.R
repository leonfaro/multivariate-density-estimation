context("MAF integration pipeline")

test_that("main pipeline runs with MAF", {
  config <- list(n_train = 100, n_test = 100, dim = 4)
  old_wd <- setwd("../..")
  source("main.R", local = TRUE)
  main()
  setwd(old_wd)
  expect_true("maf" %in% rownames(results_table))
  expect_true(all(is.finite(results_table["maf", ])))
})

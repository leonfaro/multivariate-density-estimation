test_that("main pipeline runs with MAF", {
  config <- list(n_train = 100, n_test = 100, dim = 4)
  source("main.R", local = TRUE)
  results_table <- main()
  expect_true("maf" %in% rownames(results_table))
  expect_true(all(is.finite(results_table["maf", ])))
})


library(testthat)
source("../../00_setup.R", chdir = TRUE)

test_that("all config distributions in logp table", {
  names <- sapply(config, `[[`, "distr")
  expect_true(all(names %in% names(q_supports_logp)))
})

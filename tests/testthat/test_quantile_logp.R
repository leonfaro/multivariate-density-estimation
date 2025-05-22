library(testthat)
config_choice <- 4
source("../../00_setup.R", chdir = TRUE)
source("../../01_transport_utils.R", chdir = TRUE)

test_that("logp capability table complete", {
  dnames <- vapply(config, `[[`, character(1), "distr")
  expect_true(all(!is.na(q_supports_logp[dnames])))
})

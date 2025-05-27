library(testthat)
source("../../03_baseline.R", chdir = TRUE)

config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = d$X1))
)

parse_param_spec(config, dist_registry)

test_that("param names unique", {
  expect_true(length(unique(param_names_global)) == length(param_names_global))
})

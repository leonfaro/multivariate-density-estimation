library(testthat)
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "gamma", parm  = function(d) list(shape = softplus(d$X2), rate = softplus(d$X1)))
)
source("../../00_setup.R", chdir = TRUE)

test_that("all config distributions in logp table", {
  names <- sapply(config, `[[`, "distr")
  expect_true(all(names %in% names(q_supports_logp)))
})

library(testthat)
source("../../01_data_utils.R", chdir = TRUE)

set.seed(1)
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = abs(d$X1) + 1)),
  list(distr = "beta", parm = function(d) list(shape1 = abs(d$X1) + 2,
                                               shape2 = abs(d$X1) + 3))
)

X <- SimulateData(config, 100, 123)
res <- Standardise(X, 0.8)

# Basic properties

test_that("standardisation yields zero mean and unit sd", {
  stats <- res$stats
  for (k in seq_along(stats)) {
    expect_lt(abs(mean(res$train[, k])), 1e-8)
    expect_lt(abs(sd(res$train[, k]) - 1), 1e-8)
  }
})

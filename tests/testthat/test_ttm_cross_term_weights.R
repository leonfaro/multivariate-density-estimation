source("helper_config.R")
source("../../models/ttm_cross_term.R")

test_that("Gauss-Lobatto weights sum to one", {
  for (Q in c(2, 4, 8, 12)) {
    quad <- .gauss_lobatto_01_ct(Q, degree_t = 3)
    expect_lt(abs(sum(quad$weights) - 1), 1e-12)
  }
})


old_wd <- setwd("../..")
source("main.R")

set.seed(123)
n <- 100
res <- main()
tab <- attr(res, "tab_data")
setwd(old_wd)

test_that("logL_baseline within +/-10 for N=100", {
  expect_true(all(abs(tab$true[seq_along(config)]) <= 10))
  expect_true(all(abs(tab$true_perm[seq_along(config)]) <= 10))
})

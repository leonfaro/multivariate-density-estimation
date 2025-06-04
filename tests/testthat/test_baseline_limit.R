old_wd <- setwd("../..")
source("05_main.R")

set.seed(123)
N <- 100
full_res <- main()
setwd(old_wd)

test_that("logL_baseline within +/-10 for N=100", {
  expect_true(all(abs(full_res$normal$logL_baseline) <= 10))
  expect_true(all(abs(full_res$permutation$logL_baseline) <= 10))
})

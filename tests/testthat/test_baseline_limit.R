old_wd <- setwd("../..")
source("05_main.R")

set.seed(123)
N <- 10000
full_res <- main()
setwd(old_wd)

test_that("logL_baseline within +/-10 for N=10000", {
  expect_true(all(abs(full_res$logL_baseline) <= 10))
})

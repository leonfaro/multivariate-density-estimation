old_wd <- setwd("../..")
source("main.R")

set.seed(123)
n <- 50
res <- main()
tab1 <- res$normal
tab2 <- res$perm
setwd(old_wd)

test_that("logL_baseline within +/-10 for N=50", {
  expect_true(all(abs(tab1$logL_baseline[seq_along(config)]) <= 10))
  expect_true(all(abs(tab2$logL_baseline[seq_along(config)]) <= 10))
})

test_that("main returns recorded scatter plot", {
  expect_s3_class(res$scatter_plot, "recordedplot")
})

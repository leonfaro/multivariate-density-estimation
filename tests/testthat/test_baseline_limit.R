old_wd <- setwd("../..")
source("main.R")

set.seed(123)
n <- 50
tab1 <- main()
setwd(old_wd)

vals <- as.numeric(sub(" ±.*", "", tab1$true))
test_that("logL_true within +/-10 for N=50", {
  expect_true(all(abs(vals[seq_along(config)]) <= 10))
})

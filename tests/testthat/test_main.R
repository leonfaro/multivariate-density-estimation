old_wd <- setwd("../..")
source("05_main.R")

set.seed(42)
G <- list(N = 20, config = config, seed = 42, split_ratio = 0.7)

expect_table <- data.frame(
  dim = seq_along(G$config),
  distribution = sapply(G$config, `[[`, "distr"),
  logL_baseline = NA_real_
)

# run main with reduced N
N <- 20
res <- main()
setwd(old_wd)

test_that("main outputs baseline table", {
  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("dim", "distribution", "logL_baseline"))
  expect_equal(res$dim, expect_table$dim)
  expect_equal(res$distribution, expect_table$distribution)
  expect_true(all(is.finite(res$logL_baseline)))
})

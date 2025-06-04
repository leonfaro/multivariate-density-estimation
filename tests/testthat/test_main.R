old_wd <- setwd("../..")
source("05_main.R")

set.seed(42)
G <- list(N = 20, config = config, seed = 42, split_ratio = 0.5)

expect_table <- data.frame(
  dim = c(as.character(seq_along(G$config)), "k"),
  distribution = c(sapply(G$config, `[[`, "distr"), NA_character_),
  logL_baseline = NA_real_,
  logL_trtf = NA_real_,
  logL_ks = NA_real_
)

# run main with reduced N
N <- 20
res <- main()
setwd(old_wd)

test_that("main outputs baseline table", {
  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("dim", "distribution", "logL_baseline", "logL_trtf", "logL_ks"))
  expect_equal(res$dim, expect_table$dim)
  expect_equal(res$distribution, expect_table$distribution)
  expect_true(all(is.finite(res$logL_baseline[1:length(G$config)])))
  expect_equal(
    res$logL_baseline[length(G$config) + 1],
    sum(res$logL_baseline[1:length(G$config)])
  )
  expect_equal(
    res$logL_trtf[length(G$config) + 1],
    sum(res$logL_trtf[1:length(G$config)])
  )
  expect_equal(
    res$logL_ks[length(G$config) + 1],
    sum(res$logL_ks[1:length(G$config)])
  )
})

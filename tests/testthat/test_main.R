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

test_that("main outputs tables for both orders", {
  expect_type(res, "list")
  expect_named(res, c("normal", "permutation"))

  for (tab in res) {
    expect_s3_class(tab, "data.frame")
    expect_equal(colnames(tab), c("dim", "distribution", "logL_baseline", "logL_trtf", "logL_ks"))
    expect_equal(tab$dim, expect_table$dim)
    expect_equal(tab$distribution, expect_table$distribution)
    expect_true(all(is.finite(tab$logL_baseline[1:length(G$config)])))
    expect_equal(
      tab$logL_baseline[length(G$config) + 1],
      sum(tab$logL_baseline[1:length(G$config)])
    )
    expect_equal(
      tab$logL_trtf[length(G$config) + 1],
      sum(tab$logL_trtf[1:length(G$config)])
    )
    expect_equal(
      tab$logL_ks[length(G$config) + 1],
      sum(tab$logL_ks[1:length(G$config)])
    )
  }
})

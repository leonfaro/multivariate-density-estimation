old_wd <- setwd("../..")
source("main.R")

set.seed(42)
G <- list(n = 20, config = config, seed = 42, split_ratio = 0.5)

expect_table <- data.frame(
  dim = c(as.character(seq_along(G$config)), "k"),
  distribution = c(sapply(G$config, `[[`, "distr"), NA_character_),
  logL_baseline = NA_real_,
  logL_trtf = NA_real_,
  logL_ks = NA_real_
)

# run main with reduced n
n <- 20
res <- main()
setwd(old_wd)

test_that("main outputs combined kable table", {
  expect_s3_class(res, "knitr_kable")
  tab_data <- attr(res, "tab_data")
  expect_s3_class(tab_data, "data.frame")
  expect_true(all(is.finite(tab_data$true[1:length(G$config)])))
})

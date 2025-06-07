old_wd <- setwd("../..")
source("EDA.R")
source("00_globals.R")
source("01_data_generation.R")

set.seed(1)
G <- setup_global()
G$n <- 10
X <- gen_samples(G)


test_that("create_EDA_report produces pdf", {
  pdf_file <- tempfile(fileext = ".pdf")
  res <- create_EDA_report(X, G$config, pdf_file)
  expect_true(file.exists(res))
  unlink(res)
})

test_that("create_EDA_report prints table", {
  pdf_file <- tempfile(fileext = ".pdf")
  df <- data.frame(dim = "1", distr = "gauss", val = 0.5)
  kbl_obj <- knitr::kable(df)
  attr(kbl_obj, "tab_data") <- df
  out <- capture.output(res <- create_EDA_report(X, G$config, pdf_file,
                                                table_kbl = kbl_obj))
  expect_true(any(grepl("dim", out)))
  expect_true(file.exists(res))
  unlink(res)
})

setwd(old_wd)

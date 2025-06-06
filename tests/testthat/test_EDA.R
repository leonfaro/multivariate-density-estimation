old_wd <- setwd("../..")
source("EDA.R")
source("00_globals.R")
source("01_data_generation.R")

set.seed(1)
G <- setup_global()
G$N <- 10
X <- gen_samples(G)


test_that("create_EDA_report produces pdf", {
  pdf_file <- tempfile(fileext = ".pdf")
  res <- create_EDA_report(X, G$config, pdf_file)
  expect_true(file.exists(res))
  unlink(res)
})

setwd(old_wd)

old_wd <- setwd("../..")
source("EDA.R")
source("00_globals.R")
source("01_data_generation.R")

set.seed(1)
G <- setup_global()
G$n <- 10
X <- gen_samples(G)


test_that("create_EDA_report liefert Plot-Objekte", {
  scat <- list(
    ld_base = rnorm(5),
    ld_trtf = rnorm(5),
    ld_ks = rnorm(5),
    ld_base_p = rnorm(5),
    ld_trtf_p = rnorm(5),
    ld_ks_p = rnorm(5)
  )
  res <- create_EDA_report(X, G$config, scatter_data = scat)
  expect_type(res, "list")
  expect_true(is.list(res$plots))
  expect_s3_class(res$plots[[1]], "ggplot")
})

test_that("create_EDA_report gibt Tabelle zurueck", {
  df <- data.frame(dim = "1", distr = "gauss", val = 0.5)
  kbl_obj <- knitr::kable(df)
  attr(kbl_obj, "tab_data") <- df
  res <- create_EDA_report(X, G$config, table_kbl = kbl_obj)
  expect_s3_class(res$table, "knitr_kable")
  expect_identical(attr(res$table, "tab_data"), df)
})

test_that("create_EDA_report erstellt Parameterhistogramme", {
  Xdat <- gen_samples(G, return_params = TRUE)
  res <- create_EDA_report(Xdat$X, G$config, param_list = Xdat$params)
  expect_true(length(res$param_plots) > 0)
  expect_s3_class(res$param_plots[[1]], "ggplot")
})

setwd(old_wd)

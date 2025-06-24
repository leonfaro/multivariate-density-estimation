old_wd <- setwd("../..")
source("EDA.R")
source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ks_model.R")

set.seed(1)
perm <- c(3,4,1,2)

pd <- prepare_data(20, config, perm)
models <- fit_models(pd, config)

context("pipeline helpers")

test_that("prepare_data erstellt korrekte Listen", {
  expect_true(is.list(pd))
  expect_true(is.matrix(pd$X_tr))
  expect_equal(ncol(pd$X_tr_p), length(perm))
  expect_identical(pd$X_te_p[,1], pd$X_te[,perm[1]])
  expect_true(is.list(pd$param_list))
})

test_that("fit_models liefert Modelle und Laufzeiten", {
  expect_true(is.list(models$normal))
  expect_true(all(c("true","trtf","ks") %in% names(models$normal)))
  expect_true(all(models$runtime_normal >= 0))
})

tbls <- calc_loglik_tables(models, pd)

test_that("calc_loglik_tables gibt Tabellen zurueck", {
  expect_s3_class(tbls$tab_normal, "data.frame")
  expect_equal(nrow(tbls$tab_normal), length(config) + 1)
  expect_true(all(is.finite(tbls$tab_normal$logL_trtf)))
})

scat <- make_scatter_data(models, pd)

test_that("make_scatter_data produziert numerische Vektoren", {
  expect_true(is.numeric(scat$ld_trtf))
  expect_length(scat$ld_trtf, nrow(pd$X_te))
  expect_true(!anyNA(scat$ld_ks_p))
})

setwd(old_wd)

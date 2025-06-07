source("../../00_globals.R")
source("../../01_data_generation.R")
source("../../02_split.R")

set.seed(100)
G <- setup_global()
X <- gen_samples(G)
S <- train_test_split(X, G$split_ratio, G$seed)

N_tr <- floor(G$split_ratio * G$n)
N_te <- G$n - N_tr

context("data_split")

test_that("train_test_split respects ratio and seed", {
  expect_equal(nrow(S$X_tr), N_tr)
  expect_equal(nrow(S$X_te), N_te)
  expect_equal(sort(c(S$X_tr[,1], S$X_te[,1])), sort(X[,1]))
})


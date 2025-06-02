source("../../00_globals.R")
source("../../01_data_generation.R")
source("../../models/true_model.R")
source("../../04_evaluation.R")

set.seed(123)
G <- setup_global()
G$N <- 50
X <- gen_samples(G)
N_tr <- floor(G$split_ratio * G$N)
X_tr <- X[seq_len(N_tr), ]
X_te <- X[(N_tr + 1):nrow(X), ]

M_TRUE <- fit_TRUE(X_tr, X_te, G$config)

models <- setNames(list(M_TRUE), "TRUE")


test_that("evaluate_all returns sorted table", {
  P <- evaluate_all(X_te, models)
  expect_s3_class(P, "data.frame")
  expect_equal(colnames(P), c("model_id", "neg_logL"))
  expect_equal(nrow(P), length(models))
  expect_true(all(is.finite(P$neg_logL)))
  expect_true(all(diff(P$neg_logL) >= 0))
})

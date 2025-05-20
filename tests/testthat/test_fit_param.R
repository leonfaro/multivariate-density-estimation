library(testthat)

set.seed(123)
Sys.setenv(N_train = 500, N_test = 500)
config_choice <- 4
source("../../00_setup.R", chdir = TRUE)
source("../../01_transport_utils.R", chdir = TRUE)
source("../../02_generate_data.R", chdir = TRUE)
source("../../03_param_baseline.R", chdir = TRUE)

res <- fit_param(X_pi_train, X_pi_test, config)
mean_delta <- mean(abs(res$ll_delta_df_test$delta_ll_param))

test_that("parametric fit delta_ll finite", {
  expect_true(is.finite(mean_delta))
})


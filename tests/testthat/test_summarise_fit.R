library(testthat)
set.seed(42)
Sys.setenv(N_train = 50, N_test = 50)
config_choice <- 3
source("../../00_setup.R", chdir = TRUE)
data <- generate_data()
res <- fit_param(data$train$sample$X_pi, data$test$sample$X_pi, config)
sumtab <- summarise_fit(res$param_est, data$test$sample$X_pi, res$ll_delta_df_test, config)

required_cols <- c("mean_param1", "mean_param2", "mle_param1", "mle_param2")

test_that("summary table structure", {
  expect_true(all(required_cols %in% names(sumtab)))
  expect_equal(nrow(sumtab), length(config))
})

idx <- which(sapply(config, `[[`, "distr") == "exp")
if (length(idx) > 0) {
  test_that("single parameter handled", {
    expect_identical(sumtab$mean_param2[idx], "none")
    expect_identical(sumtab$mle_param2[idx], "none")
  })
}


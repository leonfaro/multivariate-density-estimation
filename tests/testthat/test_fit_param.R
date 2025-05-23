library(testthat)

set.seed(123)
Sys.setenv(N_total = 1000)
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "gamma", parm  = function(d) list(shape = softplus(d$X2), rate = softplus(d$X1)))
)
source("../../00_setup.R", chdir = TRUE)
data <- generate_data()
X_pi_train <- data$train$sample$X_pi
X_pi_test  <- data$test$sample$X_pi
res <- fit_param(X_pi_train, X_pi_test, config)
mean_delta <- mean(abs(res$ll_delta_df_test$delta_ll_param_avg))

test_that("parametric fit delta_ll finite", {
  expect_true(is.finite(mean_delta))
})


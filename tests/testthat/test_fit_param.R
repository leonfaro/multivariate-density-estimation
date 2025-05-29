library(testthat)

set.seed(123)
Sys.setenv(N_total = 1000)
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = d$X1)),
  list(distr = "gamma", parm  = function(d) list(shape = d$X2, rate = d$X1))
)
source("../../00_setup.R", chdir = TRUE)
data <- generate_data()
X_pi_train <- data$train$sample$X_pi
X_pi_test  <- data$test$sample$X_pi
res <- fit_joint_param(X_pi_train, X_pi_test, config)
mean_ll <- mean(res$ll_df_test$ll_true)

test_that("parametric fit loglik finite", {
  expect_true(is.finite(mean_ll))
})


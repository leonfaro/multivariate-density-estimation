N <- 1000
Sys.setenv(N_train = N, N_test = N)

source("00_setup.R")
source("01_transport_utils.R")
source("02_generate_data.R")
source("03_param_baseline.R")

param_res <- fit_param(X_pi_train, X_pi_test, config)
param_est <- param_res$param_est
ll_delta_df_test <- param_res$ll_delta_df_test
ll_delta_df_test$delta_ll_param <-
  ll_delta_df_test$ll_true_sum - ll_delta_df_test$ll_param_sum

## mean parameter values on test data
mean_param_test <- sapply(seq_len(K), function(k) {
  if (k == 1) {
    param_est[[1]]$mean
  } else if (k == 2) {
    rate <- exp(param_est[[2]]$a + param_est[[2]]$b * X_pi_test[, 1])
    mean(rate)
  } else if (k == 3) {
    shape <- exp(param_est[[3]]$a) * X_pi_test[, 2]
    mean(shape)
  } else {
    NA_real_
  }
})

## first MLE component from parametric fit
mle_param <- sapply(seq_len(K), function(k) {
  unlist(param_est[[k]])[1]
})

ll_delta_df_test$mean_param_test <- round(mean_param_test, 3)
ll_delta_df_test$mle_param <- round(mle_param, 3)

print(ll_delta_df_test[
  , c(
    "dim", "distribution", "ll_true_sum", "ll_param_sum", "delta_ll_param",
    "mean_param_test", "mle_param"
  )
])

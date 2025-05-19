N <- 1000
Sys.setenv(N_train = N, N_test = N)

source("00_setup.R")
source("01_transport_utils.R")
source("02_generate_data.R")
source("03_param_baseline.R")

param_res <- fit_param(X_pi_train, X_pi_test, config)
ll_delta_df_test <- param_res$ll_delta_df_test
ll_delta_df_test$delta_ll_param <-
  ll_delta_df_test$ll_true_sum - ll_delta_df_test$ll_param_sum

print(ll_delta_df_test[
  , c("dim", "distribution", "ll_true_sum", "ll_param_sum", "delta_ll_param")
])

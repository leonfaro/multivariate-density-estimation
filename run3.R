N <- 100
# Choose which configuration to use: 3 or 4
config_choice <- 3
Sys.setenv(N_train = N, N_test = N)

source("00_setup.R")
source("01_transport_utils.R")
source("02_generate_data.R")
source("03_param_baseline.R")

param_res <- fit_param(X_pi_train, X_pi_test, config)
param_est <- param_res$param_est
ll_delta_df_test <- summarise_fit(param_est, X_pi_test, param_res$ll_delta_df_test)

ll_delta_df_test$delta_ll_param <-
  ll_delta_df_test$ll_true_sum - ll_delta_df_test$ll_param_sum

print(ll_delta_df_test[
  , c(
    "dim", "distribution", "ll_true_sum", "ll_param_sum", "delta_ll_param",
    "mean_param_test", "mle_param"
  )
])

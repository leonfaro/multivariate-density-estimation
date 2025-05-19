N <- 50
Sys.setenv(N_train = N, N_test = N)

source("00_setup.R")
source("01_transport_utils.R")
source("02_generate_data.R")
source("03_param_baseline.R")

param_res <- fit_param(X_pi_train, X_pi_test, config)
ll_delta_df_test <- param_res$ll_delta_df_test

print(ll_delta_df_test[, 1:4])

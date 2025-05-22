source("00_setup.R")
set.seed(SEED)
source("01_transport_utils.R")
Sys.setenv(N_train = 500, N_test = 500)
source("02_generate_data.R")
source("03_param_baseline.R")

X_train <- X_pi_train[, 1:2]
X_test  <- X_pi_test[, 1:2]
config_bi <- config[1:2]

res <- fit_param(X_train, X_test, config_bi)
summary_tab <- summarise_fit(res$param_est, X_test, res$ll_delta_df_test, config_bi)
print(summary_tab)

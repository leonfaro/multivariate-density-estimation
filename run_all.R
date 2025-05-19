source("00_setup.R")
source("01_transport_utils.R")
source("02_generate_data.R")
source("03_param_baseline.R")
param_res <- fit_param(X_pi_train, X_pi_test, config)
param_est <- param_res$param_est
ll_delta_df_test <- param_res$ll_delta_df_test
true_ll_mat_test <- param_res$true_ll_mat_test
param_ll_mat_test <- param_res$param_ll_mat_test

source("04_forest_models.R")
forest_res <- fit_forest(X_pi_train, X_pi_test)
model  <- forest_res$model
LD_hat <- forest_res$LD_hat
source("06_kernel_smoothing.R")
source("07_dvine_copula.R")
source("05_joint_evaluation.R")

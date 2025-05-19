N <- 50
Sys.setenv(N_train = N, N_test = N)

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

source("05_joint_evaluation.R")

# summarize mismatches after all computations
gauss_sum <- sum(dnorm(Z_eta_test, log = TRUE))
target_ll <- sum(ll_test) - gauss_sum
forest_mismatch <- sum(loglik_forest) - target_ll
kernel_mismatch <- sum(loglik_kernel) - target_ll
copula_mismatch <- ll_dvine_sum[length(ll_dvine_sum)] - target_ll

cat("forest logL mismatch =", round(forest_mismatch, 3), "\n")
cat("Kernel-smooth logL mismatch =", round(kernel_mismatch, 3), "\n")
cat("Copula log-likelihood mismatch =", round(copula_mismatch, 3), "\n")

eval_tab <- read.csv("results/evaluation_summary.csv")
print(eval_tab)

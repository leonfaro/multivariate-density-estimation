N <- 50
SEED <- 2023

Sys.setenv(N_train = N, N_test = N)

source("00_setup.R")
set.seed(SEED)  # override default seed from 00_setup.R
config <- config4
K <- length(config)

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

source("05_joint_evaluation.R")

# summarize mismatches after all computations
target_ll <- sum(ll_test)
forest_mismatch <- sum(loglik_trtf) - target_ll
kernel_mismatch <- sum(loglik_kernel) - target_ll

cat("trtf logL mismatch =", round(forest_mismatch, 3), "\n")
cat("kernel logL mismatch =", round(kernel_mismatch, 3), "\n")

eval_tab <- read.csv("results/evaluation_summary.csv")
print(eval_tab)


png("results/run5_logdensity_scatterplot.png")
plot(ld_hat, ld_true,
     xlab = "estimated log-density",
     ylab = "true log-density")
abline(a = 0, b = 1)
dev.off()
message("Saved plot to results/run5_logdensity_scatterplot.png")

png("results/run5_kernel_logdensity_scatterplot.png")
plot(ld_hat_ks, ld_true,
     xlab = "estimated log-density",
     ylab = "true log-density")
abline(a = 0, b = 1)
dev.off()
message("Saved plot to results/run5_kernel_logdensity_scatterplot.png")

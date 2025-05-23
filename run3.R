
N <- 50
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] == "big") N <- 10000

config_choice <- 3
Sys.setenv(N_train = N, N_test = N)

source("00_setup.R")
source("04_evaluation_diagnostics.R")
set.seed(SEED)
data <- generate_data()
write_data(data)

train_df <- data$train$df
test_df  <- data$test$df
X_pi_train <- data$train$sample$X_pi
X_pi_test  <- data$test$sample$X_pi

X_train <- as.matrix(train_df[paste0("Xpi", seq_len(K))])
X_test  <- as.matrix(test_df[paste0("Xpi", seq_len(K))])

cat("Train EDA:\n")
cat("- X_pi_train Mean: ", paste(round(colMeans(X_train), 3), collapse = ", "), "\n")
cat("- X_pi_train SD:   ", paste(round(apply(X_train, 2, sd), 3), collapse = ", "), "\n")
cat("- log_det(J)_train Range: [",
    round(min(train_df$det_J), 3), ",",
    round(max(train_df$det_J), 3), "]\n\n")

cat("Test EDA:\n")
cat("- X_pi_test Mean: ", paste(round(colMeans(X_test), 3), collapse = ", "), "\n")
cat("- X_pi_test SD:   ", paste(round(apply(X_test, 2, sd), 3), collapse = ", "), "\n")
cat("- log_det(J)_test Range: [",
    round(min(test_df$det_J), 3), ",",
    round(max(test_df$det_J), 3), "]\n\n")

summary_stats <- data.frame(
  dim = seq_len(K),
  mean_train = round(colMeans(X_train), 3),
  sd_train = round(apply(X_train, 2, sd), 3),
  mean_test = round(colMeans(X_test), 3),
  sd_test = round(apply(X_test, 2, sd), 3),
  row.names = NULL
)
print(summary_stats)

param_res <- fit_param(X_pi_train, X_pi_test, config)
param_est <- param_res$param_est
ll_delta_df_test <- summarise_fit(param_est, X_pi_test, param_res$ll_delta_df_test, config)

print(ll_delta_df_test[
  , c(
    "dim", "distribution", "ll_true_avg", "ll_param_avg", "delta_ll_param_avg",
    "mean_param_test", "mle_param"
  )
])

diagnostics_dir <- "diagnostics_output"
run_all_diagnostics(X_pi_train, X_pi_test, param_est,
                    config, dist_registry, diagnostics_dir)

source("dump_run3_code.R")

N <- 500
# fixed three-dimensional setup
config_choice <- 3
Sys.setenv(N_train = N, N_test = N)

source("00_setup.R")
source("01_transport_utils.R")
source("02_generate_data.R")

# compute exploratory summaries from the stored CSVs
train_df <- read.csv("results/train_data.csv")
test_df  <- read.csv("results/test_data.csv")

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

source("03_param_baseline.R")

param_res <- fit_param(X_pi_train, X_pi_test, config)
param_est <- param_res$param_est
ll_delta_df_test <- summarise_fit(param_est, X_pi_test, param_res$ll_delta_df_test, config)

ll_delta_df_test$delta_ll_param <-
  ll_delta_df_test$ll_true_sum - ll_delta_df_test$ll_param_sum

print(ll_delta_df_test[
  , c(
    "dim", "distribution", "ll_true_sum", "ll_param_sum", "delta_ll_param",
    "mean_param_test", "mle_param"
  )
])

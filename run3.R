
N <- 10000
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] == "big") N <- 10000

config_choice <- 3
Sys.setenv(N_train = N, N_test = N)

source("00_setup.R")
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
ll_delta_df_test <- summarise_fit(param_est, X_pi_test,
                                  param_res$ll_delta_df_test, config)

print(ll_delta_df_test[
  , c(
    "dim", "distribution", "ll_true_avg", "ll_param_avg",
    "delta_ll_param_avg", "mean_param1", "mean_param2",
    "mle_param1", "mle_param2"
  )
])

diagnostics_dir <- "diagnostics_output"
run_all_diagnostics(X_pi_train, X_pi_test, param_est,
                    config, dist_registry, diagnostics_dir)

# Basic EDA on diagnostic CSVs
diag_folder_path <- diagnostics_dir
if (K >= 3 && config[[3]]$distr == "gamma") {
  betas_csv_path <- file.path(diag_folder_path,
                              "dim3_gamma_estimated_betas.csv")
  if (file.exists(betas_csv_path)) {
    betas_table <- read.csv(betas_csv_path)
    cat("\n--- Dim 3 (Gamma) Estimated Beta Coefficients ---\n")
    print(betas_table)
  }

  params_train_csv_path <- file.path(diag_folder_path,
                                     "dim3_gamma_params_logpdf_train.csv")
  if (file.exists(params_train_csv_path)) {
    params_train_table <- read.csv(params_train_csv_path)
    cat("\n--- Dim 3 (Gamma) Diagnostics (Train Data) ---\n")
    if (all(c("true_shape", "fitted_shape") %in%
            colnames(params_train_table))) {
      mse_shape <- mean((params_train_table$true_shape -
                         params_train_table$fitted_shape)^2,
                        na.rm = TRUE)
      cat(sprintf("MSE (true_shape vs fitted_shape): %.6f\n", mse_shape))
      cat("Summary true_shape:\n")
      print(summary(params_train_table$true_shape))
      cat("Summary fitted_shape:\n")
      print(summary(params_train_table$fitted_shape))
    }
    if (all(c("true_rate", "fitted_rate") %in%
            colnames(params_train_table))) {
      mse_rate <- mean((params_train_table$true_rate -
                        params_train_table$fitted_rate)^2,
                       na.rm = TRUE)
      cat(sprintf("MSE (true_rate vs fitted_rate): %.6f\n", mse_rate))
      cat("Summary true_rate:\n")
      print(summary(params_train_table$true_rate))
      cat("Summary fitted_rate:\n")
      print(summary(params_train_table$fitted_rate))
    }
    sum_true_ll <- sum(params_train_table$true_log_pdf, na.rm = TRUE)
    sum_fitted_ll <- sum(params_train_table$fitted_log_pdf, na.rm = TRUE)
    cat(sprintf("Sum of true_log_pdf (train): %.6f\n", sum_true_ll))
    cat(sprintf("Sum of fitted_log_pdf (train): %.6f\n", sum_fitted_ll))
    delta_ll_dim3_train <- sum_true_ll - sum_fitted_ll
    cat(sprintf("Delta LL (true - fitted, train, dim3): %.6f\n",
                delta_ll_dim3_train))
  }
} else {
  message("Skipping detailed Gamma Dim 3 EDA.")
}

source("dump_run3_code.R")


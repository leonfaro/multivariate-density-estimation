source("00_setup.R")
set.seed(SEED)
Sys.setenv(N_train = 500, N_test = 500)
data <- generate_data()

X_train <- data$train$sample$X_pi[, 1:2]
X_test  <- data$test$sample$X_pi[, 1:2]
config_bi <- config[1:2]

res <- fit_param(X_train, X_test, config_bi)
summary_tab <- summarise_fit(res$param_est, X_test, res$ll_delta_df_test, config_bi)
print(summary_tab)

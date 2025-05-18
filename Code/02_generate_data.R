source("01_transport_utils.R")

set.seed(2044)
N_train <- as.integer(Sys.getenv("N_train", "500"))
samp_train <- pi_sample(N_train)
X_pi_train <- samp_train$X_pi
U_eta_train <- samp_train$U_eta
Z_eta_train <- samp_train$Z_eta
logd_train <- samp_train$logd

detJ_train <- det_J(logd_train)
ll_train <- loglik(Z_eta_train, detJ_train)
cat("Train EDA:\n")
cat("- X_pi_train Mean: ", paste(round(colMeans(X_pi_train), 3), collapse = ", "), "\n")
cat("- X_pi_train SD:   ", paste(round(apply(X_pi_train, 2, sd), 3), collapse = ", "), "\n")
cat("- det(J)_train Range: [", round(min(detJ_train), 3), ",", round(max(detJ_train), 3), "]\n\n")

set.seed(2045)
N_test <- as.integer(Sys.getenv("N_test", "500"))
samp_test <- pi_sample(N_test)
X_pi_test <- samp_test$X_pi
U_eta_test <- samp_test$U_eta
Z_eta_test <- samp_test$Z_eta
logd_test <- samp_test$logd

detJ_test <- det_J(logd_test)
ll_test <- loglik(Z_eta_test, detJ_test)
cat("Test EDA:\n")
cat("- X_pi_test Mean: ", paste(round(colMeans(X_pi_test), 3), collapse = ", "), "\n")
cat("- X_pi_test SD:   ", paste(round(apply(X_pi_test, 2, sd), 3), collapse = ", "), "\n")
cat("- det(J)_test Range: [", round(min(detJ_test), 3), ",", round(max(detJ_test), 3), "]\n\n")

if (!dir.exists("results")) dir.create("results")

# combine training data into one CSV
train_df <- data.frame(
  X_pi_train, U_eta_train, Z_eta_train, logd_train,
  check.names = FALSE
)
colnames(train_df) <- c(
  paste0("Xpi",  seq_len(K)),
  paste0("Ueta", seq_len(K)),
  paste0("Zeta", seq_len(K)),
  paste0("logd", seq_len(K))
)
write.csv(train_df, "results/train_data.csv", row.names = FALSE)

# combine test data into one CSV
test_df <- data.frame(
  X_pi_test, U_eta_test, Z_eta_test, logd_test,
  check.names = FALSE
)
colnames(test_df) <- colnames(train_df)
write.csv(test_df, "results/test_data.csv", row.names = FALSE)

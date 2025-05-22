out_file <- file.path("results", "code_and_output.txt")
writeLines(code_text, out_file)

output_lines <- c("###start_output_run3.R###", capture.output({
  cat("Train EDA:\n")
  cat("- X_pi_train Mean: ", paste(round(colMeans(X_train), 3), collapse = ", "), "\n")
  cat("- X_pi_train SD:   ", paste(round(apply(X_train, 2, sd), 3), collapse = ", "), "\n")
  cat("- log_det(J)_train Range: [", round(min(train_df$det_J), 3), ",", round(max(train_df$det_J), 3), "]\n\n")

  cat("Test EDA:\n")
  cat("- X_pi_test Mean: ", paste(round(colMeans(X_test), 3), collapse = ", "), "\n")
  cat("- X_pi_test SD:   ", paste(round(apply(X_test, 2, sd), 3), collapse = ", "), "\n")
  cat("- log_det(J)_test Range: [", round(min(test_df$det_J), 3), ",", round(max(test_df$det_J), 3), "]\n\n")

  print(summary_stats)
  print(ll_delta_df_test[
    , c("dim", "distribution", "ll_true_avg", "ll_param_avg", "delta_ll_param_avg", "mean_param_test", "mle_param")
  ])
}), "###end_output_run3.R###")

con <- file(out_file, open = "a")
writeLines(output_lines, con)
close(con)

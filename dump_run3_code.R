files <- c(
  "00_setup.R",
  "01_map_definition_S.R",
  "02_sampling.R",
  "03_baseline.R"
)
code_lines <- lapply(files, function(f) {
  lines <- readLines(f)
  lines <- gsub("#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  c(paste0("###start_", f, "###"), lines, paste0("###end_", f, "###"))
})
# extract only N assignment lines from run3.R
run3_lines <- readLines("run3.R")
run3_N_lines <- grep("N\\s*<-", run3_lines, value = TRUE)
run3_N_lines <- gsub("#.*$", "", run3_N_lines)
run3_N_lines <- trimws(run3_N_lines)
run3_N_lines <- run3_N_lines[nzchar(run3_N_lines)]
n_section <- c("###start_N_run3.R###", run3_N_lines, "###end_N_run3.R###")

code_text <- unlist(c(code_lines, list(n_section)))
if (!dir.exists("results")) dir.create("results")
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
  print(tbl_out)

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
      if (all(c("true_shape", "fitted_shape") %in% colnames(params_train_table))) {
        mse_shape <- mean((params_train_table$true_shape -
                           params_train_table$fitted_shape)^2,
                          na.rm = TRUE)
        cat(sprintf("MSE (true_shape vs fitted_shape): %.6f\n", mse_shape))
        cat("Summary true_shape:\n")
        print(summary(params_train_table$true_shape))
        cat("Summary fitted_shape:\n")
        print(summary(params_train_table$fitted_shape))
      }
      if (all(c("true_rate", "fitted_rate") %in% colnames(params_train_table))) {
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
}), "###end_output_run3.R###")

con <- file(out_file, open = "a")
writeLines(output_lines, con)
close(con)

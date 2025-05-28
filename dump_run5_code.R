files <- c(
  "00_setup.R",
  "01_map_definition_S.R",
  "02_sampling.R",
  "03_baseline.R",
  "04_forest_models.R",
  "05_joint_evaluation.R",
  "06_kernel_smoothing.R"
)
code_lines <- lapply(files, function(f) {
  lines <- readLines(f)
  lines <- gsub("#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  c(paste0("###start_", f, "###"), lines, paste0("###end_", f, "###"))
})
run5_lines <- readLines("run5.R")
run5_N_lines <- grep("N\\s*<-", run5_lines, value = TRUE)
run5_N_lines <- gsub("#.*$", "", run5_N_lines)
run5_N_lines <- trimws(run5_N_lines)
run5_N_lines <- run5_N_lines[nzchar(run5_N_lines)]
n_section <- c("###start_N_run5.R###", run5_N_lines, "###end_N_run5.R###")

code_text <- unlist(c(code_lines, list(n_section)))
if (!dir.exists("results")) dir.create("results")
out_file <- file.path("results", "code_and_output_run5.txt")
writeLines(code_text, out_file)

output_lines <- c("###start_output_run5.R###", capture.output({
  cat("Train EDA:\n")
  cat("- X_pi_train Mean: ", paste(round(colMeans(X_pi_train), 3), collapse = ", "), "\n")
  cat("- X_pi_train SD:   ", paste(round(apply(X_pi_train, 2, sd), 3), collapse = ", "), "\n")
  cat("- log_det(J)_train Range: [", round(min(train_df$det_J), 3), ",", round(max(train_df$det_J), 3), "]\n\n")

  cat("Test EDA:\n")
  cat("- X_pi_test Mean: ", paste(round(colMeans(X_pi_test), 3), collapse = ", "), "\n")
  cat("- X_pi_test SD:   ", paste(round(apply(X_pi_test, 2, sd), 3), collapse = ", "), "\n")
  cat("- log_det(J)_test Range: [", round(min(test_df$det_J), 3), ",", round(max(test_df$det_J), 3), "]\n\n")

  print(eval_tab)
}), "###end_output_run5.R###")

con <- file(out_file, open = "a")
writeLines(output_lines, con)
close(con)

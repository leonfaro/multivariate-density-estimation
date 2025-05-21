files <- c(
  "00_setup.R",
  "01_transport_utils.R",
  "02_generate_data.R",
  "03_param_baseline.R",
  "04_forest_models.R",
  "06_kernel_smoothing.R",
  "05_joint_evaluation.R",
  "run5.R"
)
code_lines <- lapply(files, function(f) {
  lines <- readLines(f)
  lines <- gsub("#.*$", "", lines)
  lines <- trimws(lines)
  lines[nzchar(lines)]
})
code_text <- unlist(code_lines)
if (!dir.exists("results")) dir.create("results")
out_file <- file.path("results", "code_and_output_run5.txt")
writeLines(code_text, out_file)

output_lines <- capture.output({
  cat("trtf logL mismatch =", round(forest_mismatch, 3), "\n")
  cat("kernel logL mismatch =", round(kernel_mismatch, 3), "\n")
  print(eval_tab)
})
con <- file(out_file, open = "a")
writeLines(output_lines, con)
close(con)

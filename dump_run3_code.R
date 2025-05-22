files <- c(
  "00_setup.R",
  "01_transport_utils.R",
  "02_generate_data.R",
  "03_param_baseline.R",
  "run3.R"
)
code_lines <- lapply(files, function(f) {
  lines <- readLines(f)
  lines <- gsub("#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  c(paste0("###start_", f, "###"), lines, paste0("###end_", f, "###"))
})
code_text <- unlist(code_lines)
if (!dir.exists("results")) dir.create("results")
out_file <- file.path("results", "code_and_output.txt")
writeLines(code_text, out_file)

output_lines <- c("###start_output_run3.R###",
  capture.output(source("run3.R", echo = FALSE)),
  "###end_output_run3.R###")
con <- file(out_file, open = "a")
writeLines(output_lines, con)
close(con)

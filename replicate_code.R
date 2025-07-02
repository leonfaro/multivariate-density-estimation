extract_sources <- function(main_file = "main.R") {
  lines <- readLines(main_file, warn = FALSE)
  pattern <- "source\(\"([^"]+)\"\)"
  matches <- regmatches(lines, gregexpr(pattern, lines, perl = TRUE))
  src_files <- unlist(lapply(matches, function(x) if (length(x) > 0) gsub(pattern, "\\1", x)))
  unique(src_files)
}

replicate_code_scripts <- function(main_file = "main.R", outfile = "replicated_code.txt") {
  src_files <- extract_sources(main_file)
  output_lines <- character()
  for (f in src_files) {
    if (file.exists(f)) {
      lines <- readLines(f, warn = FALSE)
      lines <- sub("#.*$", "", lines)
      lines <- trimws(lines)
      lines <- lines[nchar(lines) > 0]
      output_lines <- c(output_lines, paste0("### Begin ", f, " ###"), lines, paste0("### End ", f, " ###"), "")
    }
  }
  writeLines(output_lines, outfile)
}

if (sys.nframe() == 0L) {
  replicate_code_scripts()
}

#!/usr/bin/env Rscript

# Export all TTM R sources into a single text file with per-file headers.
# Usage examples:
#   Rscript scripts/export_ttm_sources.R
#   Rscript scripts/export_ttm_sources.R --out=ttm_all_scripts.txt
#   Rscript scripts/export_ttm_sources.R --src=models/ttm --out=out/ttm_sources.txt

args <- commandArgs(trailingOnly = TRUE)

print_help <- function() {
  cat("Export TTM sources to a single text file\n")
  cat("\nOptions:\n")
  cat("  --src=DIR   Source directory (default: models/ttm)\n")
  cat("  --out=FILE  Output file (default: ttm_all_scripts.txt in repo root)\n")
  cat("  -h, --help  Show this help\n")
}

if (length(args) && any(args %in% c("-h", "--help"))) {
  print_help(); quit(status = 0)
}

# Try to find the repository root (contains 00_globals.R)
get_root <- function() {
  if (exists("root_path", inherits = FALSE)) return(root_path)
  p <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (file.exists(file.path(p, "00_globals.R"))) return(p)
    np <- dirname(p); if (identical(np, p)) break; p <- np
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

root <- get_root()

# Defaults
src_dir <- file.path(root, "models", "ttm")
out_file <- file.path(root, "ttm_all_scripts.txt")

# Parse key=value args
if (length(args)) {
  for (a in args) {
    if (grepl("^--src=", a)) src_dir <- normalizePath(sub("^--src=", "", a), winslash = "/", mustWork = FALSE)
    if (grepl("^--out=", a)) {
      of <- sub("^--out=", "", a)
      # If relative, interpret relative to root
      if (!grepl("^/|^[A-Za-z]:", of)) of <- file.path(root, of)
      out_file <- normalizePath(of, winslash = "/", mustWork = FALSE)
    }
  }
}

if (!dir.exists(src_dir)) {
  stop(sprintf("Source directory not found: %s", src_dir))
}

files <- list.files(src_dir, pattern = "\\.R$", full.names = TRUE)
files <- sort(files)
if (length(files) == 0L) stop(sprintf("No R files found in %s", src_dir))

# Ensure output directory exists
out_dir <- dirname(out_file)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

con <- file(out_file, open = "wt")
on.exit(try(close(con), silent = TRUE), add = TRUE)

for (f in files) {
  bn <- basename(f)
  writeLines(sprintf("# %s", bn), con)
  txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  if (length(txt)) writeLines(txt, con)
  writeLines(c("", ""), con) # separator
}

cat(sprintf("Wrote %d files to %s\n", length(files), out_file))


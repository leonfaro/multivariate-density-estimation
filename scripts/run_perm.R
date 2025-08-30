#!/usr/bin/env Rscript
## Helper: run main.R with custom permutation
## Usage: Rscript scripts/run_perm.R 4,1,2,3

args <- commandArgs(trailingOnly = TRUE)
perm_str <- if (length(args) >= 1L) args[[1]] else "4,1,2,3"
perm_vals <- as.integer(strsplit(perm_str, ",", fixed = TRUE)[[1]])

source("main.R")
perm <- perm_vals
tab <- main()
print(tab)
if (exists("timing_table", inherits = FALSE)) print(timing_table)

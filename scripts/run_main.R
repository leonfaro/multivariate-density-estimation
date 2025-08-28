#!/usr/bin/env Rscript
options(error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 1) })
options(mc.cores = min(parallel::detectCores(), 10L))
n_override <- as.integer(Sys.getenv("N_OVERRIDE", "100"))

source("00_globals.R")
source("main.R")
assign("n", n_override, envir = .GlobalEnv)

tab <- main()
print(tab)
cat("[INFO] Pipeline executed (train/test only)\n")

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
tryCatch({
  write.csv(tab, file = "results/last_run_table.csv", row.names = FALSE)
  saveRDS(tab, file = "results/last_run_table.rds")
}, error = function(e) message("WARN (save): ", conditionMessage(e)))

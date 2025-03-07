# main.R

setwd("/Users/leonkiafaro/Documents/multivariate-density-estimation/Code/2025_03_07")

source("00_setup.R")
source("01_distribution.R")
source("02_dgp.R")

d_used <- 25
N      <- 100
debug_mode <- FALSE

all_data <- generate_data(
  N           = N,
  d           = d_used,
  config      = config[1:d_used],
  distLibrary = distLibrary,
  debug       = debug_mode
)

cat("\nHead(all_data):\n")
print(head(all_data))

cat("\nCheck for NA values:\n")
na_ix <- which(is.na(all_data), arr.ind=TRUE)
print(na_ix)

logdens_vec <- compute_logdensity(
  Y           = all_data,
  d           = d_used,
  config      = config[1:d_used],
  distLibrary = distLibrary,
  debug       = debug_mode
)

cat("\nMean log-density =", mean(logdens_vec), "\n")



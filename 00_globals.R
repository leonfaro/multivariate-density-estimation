library(dplyr)
library(parallel)
library(tram)
library(trtf)
if (!exists("NC")) {
  NC <- suppressWarnings(as.integer(parallel::detectCores()))
  if (!is.finite(NC) || is.na(NC) || NC < 1L) NC <- 1L
}
options(mc.cores = NC)

softplus <- function(x) log1p(exp(x))


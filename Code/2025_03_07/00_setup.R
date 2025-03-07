# 00_setup.R

library(tram)       # diverse transformations
library(trtf)       # transformation forests
library(ks)         # kernel smoothing
library(partykit)   # partykit
library(party)      # party
library(foreach)    # parallel
library(doParallel) # parallel

set.seed(123)       # fix seed
options(digits=4)   # 4 decimals

init_parallel <- function(ncores=2) {
  cl <- makeCluster(ncores)    # cluster creation
  registerDoParallel(cl)       # register
  cl                           # return cluster
}





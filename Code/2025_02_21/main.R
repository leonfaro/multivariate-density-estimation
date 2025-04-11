# main.R

setwd("/Users/leonkiafaro/Documents/multivariate-density-estimation/Code/2025_02_21")

source("00_setup.R")
source("01_distribution.R")  
source("02_dgp.R")           
source("03_models.R")        
source("04_analyze.R")       


all_data <- my_dgp_d(N=50, d=3, models=models, cond=cond,
                     chunk_size=NULL, parallel=FALSE, debug=FALSE)

idx <- sample(nrow(all_data), 70)
trainD <- all_data[idx, ]
testD  <- all_data[-idx, ]

res <- analyze_models(trainD, dtest=testD)



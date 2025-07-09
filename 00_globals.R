
library(parallel)
if (!exists("NC")) NC <- detectCores()
options(mc.cores = NC)

softplus <- function(x) log1p(exp(x))



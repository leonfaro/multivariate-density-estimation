config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "gamma", parm  = function(d) list(shape = softplus(d$X2), rate = softplus(d$X1)))
)

source("00_setup.R")
set.seed(123)
dat <- generate_data(N_total = 10)
stopifnot(!any(is.na(dat$train$df)))
stopifnot(!any(is.na(dat$test$df)))
cat("basic tests passed\n")

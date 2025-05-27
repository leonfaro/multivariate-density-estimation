config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = d$X1)),
  list(distr = "beta", parm  = function(d) list(shape1 = d$X2, shape2 = 1))
)

N <- 1000
source("00_setup.R")
set.seed(SEED)
data <- pi_sample(N, config)

cdf_vals <- numeric(N)
logcdf_vals <- numeric(N)
for (i in seq_len(N)) {
  x_prev <- data$X_pi[i, 1:2]
  xk <- data$X_pi[i, 3]
  cdf_vals[i] <- cdf_k(3, xk, x_prev, config, log = FALSE)
  logcdf_vals[i] <- cdf_k(3, xk, x_prev, config, log = TRUE)
}

res <- list(
  cdf_range = range(cdf_vals),
  logcdf_range = range(logcdf_vals)
)
str(res)

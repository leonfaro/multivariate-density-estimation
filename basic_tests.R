source("00_setup.R")
source("01_transport_utils.R")

set.seed(1)
samp <- pi_sample(10)
stopifnot(dim(samp$X_pi)[1] == 10)
stopifnot(dim(samp$X_pi)[2] == K)
stopifnot(!any(is.na(samp$X_pi)))

ld <- matrix(0, nrow=5, ncol=K)
z  <- matrix(rnorm(5*K), nrow=5)
ll <- loglik(z, ld)
stopifnot(is.numeric(ll))
stopifnot(all(is.finite(ll)))

cat("All basic tests passed.\n")

library(testthat)
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = d$X1)),
  list(distr = "gamma", parm  = function(d) list(shape = d$X2, rate = d$X1))
)
source("../../00_setup.R", chdir = TRUE)

# Softplus positivity
vals <- seq(-100, 100, length.out = 201)

test_that("softplus positive", {
  expect_true(all(softplus(vals) > 0))
})

# Synthetic gamma fit
set.seed(1)
N <- 10000
x_prev <- rnorm(N)
true_theta <- c(0.5, -0.2, -0.1, 0.3)
shape <- softplus(true_theta[1] + true_theta[2] * x_prev)
rate  <- softplus(true_theta[3] + true_theta[4] * x_prev)
xs <- rgamma(N, shape = shape, rate = rate)

nll <- make_generalized_nll("gamma", matrix(x_prev, ncol = 1), xs)
fit <- safe_optim(rep(0, 4), nll)

shape_hat <- softplus(fit$par[1] + fit$par[2] * x_prev)
rate_hat  <- softplus(fit$par[3] + fit$par[4] * x_prev)
rel <- c(mean(abs(shape_hat - shape) / shape),
         mean(abs(rate_hat - rate) / rate))

test_that("gamma parameters recovered", {
  expect_true(all(rel < 0.01))
})

test_that("nll monotonicity", {
  expect_lte(nll(fit$par), nll(rep(0, 4)))
})

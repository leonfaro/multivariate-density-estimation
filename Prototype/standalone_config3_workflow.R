# Standalone workflow based on config3
# Implements target sampling, triangular map and density evaluation

# --- Configuration -------------------------------------------------------
config3 <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = d$X1)),
  list(distr = "gamma", parm  = function(d) list(shape = d$X2, rate = 1))
)

# --- Helper Functions ----------------------------------------------------
set.seed(42)
K <- 3
n_train <- 50
n_test  <- 50

lambda_fun <- function(x1) abs(x1) + 1e-8

sample_pi <- function(n) {
  X <- matrix(NA_real_, nrow = n, ncol = K,
              dimnames = list(NULL, c("X1", "X2", "X3")))
  for (i in seq_len(n)) {
    x1 <- rnorm(1)
    lam <- lambda_fun(x1)
    x2 <- rexp(1, rate = lam)
    x3 <- rgamma(1, shape = x2, rate = 1)
    X[i, ] <- c(x1, x2, x3)
  }
  X
}

eta_density <- function(z) prod(dnorm(z))

sample_eta <- function(n) {
  matrix(rnorm(n * K), nrow = n, ncol = K,
         dimnames = list(NULL, c("Z1", "Z2", "Z3")))
}

S1 <- function(x1) qnorm(pnorm(x1))
S2 <- function(x1, x2) {
  qnorm(1 - exp(-lambda_fun(x1) * x2))
}
S3 <- function(x1, x2, x3) {
  qnorm(pgamma(x3, shape = x2, rate = 1))
}

S <- function(x) {
  cbind(S1(x[,1]), S2(x[,1], x[,2]), S3(x[,1], x[,2], x[,3]))
}

S1_inv <- function(z1) z1
S2_inv <- function(z2, x1) {
  -log1p(-pnorm(z2)) / lambda_fun(x1)
}
S3_inv <- function(z3, x1, x2) {
  qgamma(pnorm(z3), shape = x2, rate = 1)
}

S_inv <- function(z) {
  n <- nrow(z)
  X <- matrix(NA_real_, nrow = n, ncol = K,
              dimnames = list(NULL, c("X1", "X2", "X3")))
  for (i in seq_len(n)) {
    x1 <- S1_inv(z[i,1])
    x2 <- S2_inv(z[i,2], x1)
    x3 <- S3_inv(z[i,3], x1, x2)
    X[i, ] <- c(x1, x2, x3)
  }
  X
}

# Density evaluation ------------------------------------------------------
pi_density <- function(x) {
  z <- S(x)
  phi_z <- dnorm(z)
  lambda <- lambda_fun(x[,1])
  J2 <- lambda * exp(-lambda * x[,2]) / dnorm(z[,2])
  J3 <- x[,3]^(x[,2]-1) * exp(-x[,3]) / (gamma(x[,2]) * dnorm(z[,3]))
  dens <- dnorm(z[,1]) * dnorm(z[,2]) * dnorm(z[,3]) * J2 * J3
  dens
}

# Conditional sampling ----------------------------------------------------
sample_conditional <- function(n, x_fixed = list(X1 = NULL, X2 = NULL, X3 = NULL)) {
  Z <- sample_eta(n)
  X <- matrix(NA_real_, nrow = n, ncol = K,
              dimnames = list(NULL, c("X1", "X2", "X3")))
  for (i in seq_len(n)) {
    x1 <- if (is.null(x_fixed$X1)) S1_inv(Z[i,1]) else x_fixed$X1
    x2 <- if (is.null(x_fixed$X2)) S2_inv(Z[i,2], x1) else x_fixed$X2
    x3 <- if (is.null(x_fixed$X3)) S3_inv(Z[i,3], x1, x2) else x_fixed$X3
    X[i, ] <- c(x1, x2, x3)
  }
  X
}

# --- Generate Training and Test Data ------------------------------------
X_train <- sample_pi(n_train)
X_test  <- sample_pi(n_test)

# Simple checks -----------------------------------------------------------
stopifnot(!any(is.na(X_train)))
stopifnot(!any(is.na(X_test)))

summary(X_train)
summary(X_test)

# Example usage: density of first training observation
pi_density(matrix(X_train[1, ], nrow = 1))


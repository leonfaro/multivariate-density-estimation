#!/usr/bin/env Rscript
# Minimal smoke test for ttm_cross_term: forward-KL with fixed B-spline time basis

options(mde.verbose = TRUE)
set.seed(123)

# Source module
root <- getwd()
source(file.path(root, "00_globals.R"))
source(file.path(root, "models/ttm_cross_term.R"))

# Generate small 2D dataset
N <- 50L; K <- 2L
X <- matrix(rnorm(N * K), ncol = K)

# Train/test split
idx <- sample.int(N)
n_tr <- floor(0.7 * N)
X_tr <- X[idx[seq_len(n_tr)], , drop = FALSE]
X_te <- X[idx[(n_tr + 1L):N], , drop = FALSE]

S <- list(X_tr = X_tr, X_te = X_te)

# Train with defaults (lean B-spline path, Q=24, lambdas fixed)
fit <- trainCrossTermMap(S)
stopifnot(is.list(fit), !is.null(fit$S))
M <- fit$S

# Predict
LD_by <- predict(M, X_te, "logdensity_by_dim")
LD_j  <- predict(M, X_te, "logdensity")

# Checks
stopifnot(is.matrix(LD_by), nrow(LD_by) == nrow(X_te), ncol(LD_by) == ncol(X_te))
stopifnot(is.numeric(LD_j), length(LD_j) == nrow(X_te))
stopifnot(all(is.finite(LD_by)), all(is.finite(LD_j)))
stopifnot(max(abs(rowSums(LD_by) - LD_j)) <= 1e-10)

# Log concise summary when verbose
if (isTRUE(getOption("mde.verbose", FALSE))) {
  df <- if (!is.null(M$coeffs) && length(M$coeffs) >= 1L && !is.null(M$coeffs[[1]]$bs_spec)) M$coeffs[[1]]$bs_spec$df else NA_integer_
  cat(sprintf("[SMOKE] N=%d, K=%d, Q=%d, df=%s\n", nrow(X_te), ncol(X_te), M$Q, as.character(df)))
}

invisible(TRUE)


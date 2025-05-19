# Sequential kernel smoothing for conditional log-densities
# - Input: training data frame `X_pi_train` with `N` rows and `K` columns
# - Output: `ks_model` object; `predict()` returns an `N_test x K` matrix of log-densities
# - Algorithm:
#   * Compute Silverman's bandwidth for each column: `bw[k] = 1.06 * sd(X[, k]) * N^(-1/5)`
#   * For a new observation `x` iterate over `k = 1, \dots, K`
#     - Accumulate weights from kernels of the previous `k-1` variables
#     - Evaluate the Gaussian kernel at `x_k` and divide by `bw[k]`
#     - Obtain the conditional density by a weighted average
#     - Store the log-density via `safe_logdens()`
#   * Repeat for all rows in `newdata`

fit_kernel <- function(data) {
  N <- nrow(data)
  K <- ncol(data)
  bw <- apply(data, 2, function(x) 1.06 * sd(x) * N^(-1/5))
  structure(list(train = data, bw = bw, K = K), class = "mykernel")
}

predict.mykernel <- function(object, newdata, type = "logdensity") {
  train <- as.matrix(object$train)
  bw <- object$bw
  K <- object$K
  test <- as.matrix(newdata)
  ntest <- nrow(test)
  ld_mat <- matrix(NA_real_, ntest, K)
  for (i in seq_len(ntest)) {
    for (k in seq_len(K)) {
      log_w <- rep(0, nrow(train))
      if (k > 1) {
        for (m in seq_len(k - 1)) {
          log_w <- log_w + dnorm((test[i, m] - train[, m]) / bw[m], log = TRUE)
        }
      }
      log_k <- dnorm((test[i, k] - train[, k]) / bw[k], log = TRUE) -
        log(bw[k])
      num <- logsumexp(log_w + log_k)
      den <- logsumexp(log_w)
      dens <- exp(num - den)
      ld_mat[i, k] <- safe_logdens(dens)
    }
  }
  ld_mat
}

ks_model <- fit_kernel(as.data.frame(X_pi_train))
KS_hat <- predict(ks_model, newdata = as.data.frame(X_pi_test), type = "logdensity")

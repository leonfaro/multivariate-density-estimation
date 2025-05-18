# Kernel smoothing for sequential conditional densities
# Input: data.frame with columns x1,x2,x3
# Gaussian kernels with Silverman's rule

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
  weights <- rep(1, nrow(train))
  for (k in seq_len(K)) {
  if (k > 1) {
  for (m in seq_len(k - 1)) {
  weights <- weights * dnorm((test[i, m] - train[, m]) / bw[m])
  }
  }
  kern_vals <- dnorm((test[i, k] - train[, k]) / bw[k]) / bw[k]
  dens <- sum(weights * kern_vals) / sum(weights)
  ld_mat[i, k] <- safe_logdens(dens)
  }
  }
  ld_mat
}

ks_model <- fit_kernel(as.data.frame(X_pi_train))
KS_hat <- predict(ks_model, newdata = as.data.frame(X_pi_test), type = "logdensity")

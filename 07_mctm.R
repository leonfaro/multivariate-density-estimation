# Minimal MCTM implementation
# loglik_fn(y, theta) = log phi(A %*% h(y)) + sum_j log h'_j(y_j)

loglik_fn <- function(y, theta) {
  y <- as.matrix(y)
  K <- ncol(y)
  n <- nrow(y)
  h_y <- matrix(NA_real_, n, K)
  log_deriv <- numeric(K)
  for (j in seq_len(K)) {
    mu <- theta$h[[j]]$mu
    sd <- theta$h[[j]]$sd
    h_y[, j] <- (y[, j] - mu) / sd
    log_deriv[j] <- -log(sd)
  }
  z <- h_y %*% t(theta$A)
  log_phi <- -0.5 * rowSums(z^2) - (K / 2) * log(2 * pi)
  ll <- log_phi + rowSums(matrix(log_deriv, n, K, byrow = TRUE))
  ll
}

fit_h <- function(y) {
  list(mu = mean(y), sd = sd(y))
}

fit_model <- function(train, init = NULL) {
  train <- as.matrix(train)
  K <- ncol(train)
  h_list <- lapply(seq_len(K), function(j) fit_h(train[, j]))
  h_train <- train
  for (j in seq_len(K)) {
    h_train[, j] <- (train[, j] - h_list[[j]]$mu) / h_list[[j]]$sd
  }
  cov_hat <- cov(h_train)
  A <- t(chol(cov_hat))
  structure(list(h = h_list, A = A), class = "mctm")
}

predict.mctm <- function(object, newdata) {
  loglik_fn(newdata, object)
}

fit_mctm <- function(train, test) {
  model <- fit_model(train)
  ll_train <- predict(model, train)
  ll_test <- predict(model, test)
  list(model = model, ll_train = ll_train, ll_test = ll_test)
}

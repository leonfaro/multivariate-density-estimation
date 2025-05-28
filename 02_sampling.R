# Input: K, config, functions pdf_k/qtf_k
# Output: eta_sample(), S_inv(), pi_sample(), sample_pi_df(), generate_data(), write_data()
# Sampling utilities

# Input: log-Jacobian matrix
# Output: vector of log-determinants
logdet_J <- function(logd) rowSums(logd)

# Input: reference samples Z_eta, log-Jacobian sums
# Output: log-likelihood under eta
loglik <- function(Z_eta, logdet_J) {
  -0.5 * rowSums(Z_eta^2) - (ncol(Z_eta)/2) * log(2 * pi) + logdet_J
}

# Input: sample size N
# Output: list with U_eta and Z_eta
eta_sample <- function(N) {
  Z_eta <- matrix(rnorm(N * K), nrow = N)
  U_eta <- pnorm(Z_eta)
  list(U_eta = U_eta, Z_eta = Z_eta)
}

# Input: matrix U_eta, cfg list, optional Z_eta
# Output: list with X_pi, U_eta, Z_eta, logd
S_inv <- function(U_eta, cfg = config, Z_eta = qnorm(U_eta)) {
  n <- nrow(U_eta)
  X_pi <- Z_eta
  logd <- Z_eta
  for (i in seq_len(n)) {
    x_prev <- numeric(0)
    for (k in seq_len(K)) {
      logu <- pnorm(Z_eta[i, k], log.p = TRUE)
      x_prev[k] <- qtf_k(k, logu, x_prev, cfg, log.p = TRUE)
      logd[i, k] <- pdf_k(k, x_prev[k], x_prev[seq_len(k - 1)],
                          cfg, log = TRUE) -
                    dnorm(Z_eta[i, k], log = TRUE)
    }
    X_pi[i, ] <- x_prev
  }
  list(X_pi = X_pi, U_eta = U_eta, Z_eta = Z_eta, logd = logd)
}

# Input: sample size N, cfg list
# Output: list with X_pi, U_eta, Z_eta, logd
pi_sample <- function(N, cfg = config) {
  ref <- eta_sample(N)
  inv <- S_inv(ref$U_eta, cfg, ref$Z_eta)
  list(X_pi = inv$X_pi, U_eta = ref$U_eta,
       Z_eta = inv$Z_eta, logd = inv$logd)
}

# Input: sample size N, cfg list
# Output: list with data frame and raw sample
sample_pi_df <- function(N, cfg = config) {
  samp <- pi_sample(N, cfg)
  logdet <- logdet_J(samp$logd)
  df <- as.data.frame(cbind(samp$X_pi, samp$U_eta, samp$Z_eta, samp$logd,
                            det_J = logdet,
                            ll_true = loglik(samp$Z_eta, logdet)))
  colnames(df) <- c(paste0("Xpi",  seq_len(K)),
                    paste0("Ueta", seq_len(K)),
                    paste0("Zeta", seq_len(K)),
                    paste0("logd", seq_len(K)),
                    "det_J", "ll_true")
  attr(df, "seed") <- SEED
  list(df = df, sample = samp)
}

# Input: total size, cfg list, split ratio
# Output: train/test list of samples and data frames
generate_data <- function(N_total = as.integer(Sys.getenv("N_total", "1000")),
                          cfg = config, split_ratio = 0.7) {
  all_dat <- sample_pi_df(N_total, cfg)
  n_train <- floor(split_ratio * N_total)
  idx_train <- seq_len(n_train)
  idx_test <- setdiff(seq_len(N_total), idx_train)

  slice_sample <- function(samp, idx) {
    lapply(samp, function(x) {
      if (is.matrix(x)) x[idx, , drop = FALSE] else x[idx]
    })
  }

  train <- list(
    df = all_dat$df[idx_train, , drop = FALSE],
    sample = slice_sample(all_dat$sample, idx_train)
  )
  test <- list(
    df = all_dat$df[idx_test, , drop = FALSE],
    sample = slice_sample(all_dat$sample, idx_test)
  )
  list(train = train, test = test)
}

# Input: generated data list, target directory
# Output: CSV files on disk
write_data <- function(data, dir = "results") {
  if (!dir.exists(dir)) dir.create(dir)
  write.csv(data$train$df, file.path(dir, "train_data.csv"), row.names = FALSE)
  write.csv(data$test$df,  file.path(dir, "test_data.csv"),  row.names = FALSE)
}

# Sampling utilities

logdet_J <- function(logd) rowSums(logd)

loglik <- function(Z_eta, logdet_J) {
  -0.5 * rowSums(Z_eta^2) - (ncol(Z_eta)/2) * log(2 * pi) + logdet_J
}

eta_sample <- function(N) {
  Z_eta <- matrix(rnorm(N * K), nrow = N)
  U_eta <- pnorm(Z_eta)
  list(U_eta = U_eta, Z_eta = Z_eta)
}

S_inv <- function(U_eta, cfg = config, Z_eta = qnorm(U_eta)) {
  n <- nrow(U_eta)
  X_pi <- matrix(NA_real_, n, K)
  logd <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    x_prev <- numeric(0)
    for (k in seq_len(K)) {
      logu <- pnorm(Z_eta[i, k], log.p = TRUE)
      x_prev[k] <- qtf_k(k, logu, x_prev, cfg, log.p = TRUE)
    }
    X_pi[i, ] <- x_prev
    for (k in seq_len(K)) {
      logd[i, k] <- pdf_k(k, X_pi[i, k], X_pi[i, seq_len(k-1)], cfg, log = TRUE) -
                    dnorm(Z_eta[i, k], log = TRUE)
    }
  }
  list(X_pi = X_pi, U_eta = U_eta, Z_eta = Z_eta, logd = logd)
}

pi_sample <- function(N, cfg = config) {
  ref <- eta_sample(N)
  inv <- S_inv(ref$U_eta, cfg, ref$Z_eta)
  list(X_pi = inv$X_pi, U_eta = ref$U_eta,
       Z_eta = inv$Z_eta, logd = inv$logd)
}

sample_pi_df <- function(N, cfg = config) {
  samp <- pi_sample(N, cfg)
  logdet <- logdet_J(samp$logd)
  ll <- loglik(samp$Z_eta, logdet)
  df <- data.frame(
    samp$X_pi, samp$U_eta, samp$Z_eta, samp$logd,
    det_J = logdet, ll_true = ll,
    check.names = FALSE
  )
  colnames(df) <- c(
    paste0("Xpi",  seq_len(K)),
    paste0("Ueta", seq_len(K)),
    paste0("Zeta", seq_len(K)),
    paste0("logd", seq_len(K)),
    "det_J", "ll_true"
  )
  attr(df, "seed") <- SEED
  list(df = df, sample = samp)
}

generate_data <- function(N_train = as.integer(Sys.getenv("N_train", "500")),
                          N_test = as.integer(Sys.getenv("N_test", "500")),
                          cfg = config) {
  train <- sample_pi_df(N_train, cfg)
  test  <- sample_pi_df(N_test, cfg)
  list(train = train, test = test)
}

write_data <- function(data, dir = "results") {
  if (!dir.exists(dir)) dir.create(dir)
  write.csv(data$train$df, file.path(dir, "train_data.csv"), row.names = FALSE)
  write.csv(data$test$df,  file.path(dir, "test_data.csv"),  row.names = FALSE)
}

# Sampling utilities

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

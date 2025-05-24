## DO NOT CHANGE 00_setup.R or 01_map_definition_S.R or 02_sampling.R or 03_baseline.R or run3.R
# Sequentielle Kernel-Glättung für Logdichten
# Eingabe: Datenmatrix `X_pi_train` (N x K)
# Ausgabe: Objekt `ks_model`; `predict()` liefert N_test x K Logdichten
# Ablauf:
#   * Bandbreite bw[k] = 1.06 * sd(X[,k]) * N^(-1/5)
#   * für jedes Test-x und k:
#       - Gewichte der früheren k-1 Variablen aufsummieren
#       - Gauß-Kernel an x_k / bw[k]
#       - gewichtete Dichte berechnen
#       - Logdichte via `safe_logdens()` speichern

if (!exists("logsumexp")) {
  logsumexp <- function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
  }
}

if (!exists("safe_logdens")) {
  safe_logdens <- function(dens, eps = EPS) {
    res <- pmax(log(eps), log(dens))
    stopifnot(all(is.finite(res)))
    res
  }
}

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

      log_k <- dnorm((test[i, k] - train[, k]) / bw[k], log = TRUE) - log(bw[k])
      numer <- logsumexp(log_w + log_k)
      denom <- logsumexp(log_w)
      dens <- exp(numer - denom)

      ld_mat[i, k] <- safe_logdens(dens)
    }
  }
  ld_mat
}

ks_model <- fit_kernel(as.data.frame(X_pi_train))
KS_hat <- predict(ks_model, newdata = as.data.frame(X_pi_test), type = "logdensity")

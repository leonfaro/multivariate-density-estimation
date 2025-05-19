# ------------------------------------------------------------------------------
# Skript: 07_dvine_copula.R
# Ziel: gemeinsamer Log-Dichte-Schätzer mittels D-Vine-Copula
#
# Eingabe:
# - X_pi_train, X_pi_test: Matrizen (N × K) mit Beobachtungen im π-Raum
# - model: fit_forest()-Objekt mit bedingten Verteilungen
# - LD_hat: Log-Dichten der Marginalen für das Testset
#
# Ausgabe:
# - loglik_dvine: Vektor mit Log-Dichte-Schätzungen für X_pi_test
#
# Algorithmus:
# - mittels predict(..., type = "distribution") U_hat_train und U_hat_test berechnen
# - dvine_structure(1:K) definieren und vinecop() auf U_hat_train anwenden
# - dvinecop() liefert log-Kopuladichten für U_hat_test
# - LD_hat mit den log-Kopuladichten addieren -> loglik_dvine
#
# Variablennamen folgen Notation.md
# ------------------------------------------------------------------------------

source("04_forest_models.R")

library(rvinecopulib)

# obtain estimated conditional CDFs from the forests
U_hat_train <- predict(model, newdata = as.data.frame(X_pi_train),
                       type = "distribution")
U_hat_test  <- predict(model, newdata = as.data.frame(X_pi_test),
                      type = "distribution")

# compute D-vine joint log-likelihood for each dimension separately
loglik_dvine_list <- lapply(seq_len(K), function(k) {
  vine_struct <- dvine_structure(1:k)
  dvine_model <- vinecop(U_hat_train[, 1:k, drop = FALSE],
                         structure = vine_struct)
  log_cop <- log(dvinecop(U_hat_test[, 1:k, drop = FALSE], dvine_model))
  rowSums(LD_hat[, 1:k, drop = FALSE]) + log_cop
})
loglik_dvine_mat <- do.call(cbind, loglik_dvine_list)
stopifnot(all(is.finite(loglik_dvine_mat)))

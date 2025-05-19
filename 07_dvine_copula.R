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

# fit a D-vine copula with fixed natural order 1:K
vine_struct <- dvine_structure(1:K)
dvine_model <- vinecop(U_hat_train, structure = vine_struct)

# log-copula densities for the test sample
log_cop_test <- dvinecop(U_hat_test, dvine_model, log = TRUE)

# joint log-density estimate via copula times marginals
loglik_dvine <- rowSums(LD_hat) + log_cop_test
stopifnot(all(is.finite(loglik_dvine)))


# Daten erzeugen ---------------------------------------------------------
# Eingabe: Umgebungsvariablen `N_train`, `N_test` (Standard 500)
# Ausgabe: CSVs unter `results/` mit Spalten `Xpi*`, `Ueta*`, `Zeta*`,
#          `logd*`, `det_J`, `ll_true` (det_J = Logdeterminante)
# Ablauf:
#   1. U_eta aus Referenz ziehen
#   2. `S_inv` anwenden -> X_pi, Z_eta, logd
#   3. logdet_J = rowSums(logd); loglik berechnen
#   4. kurze Zusammenfassung ausgeben
#   5. Daten als CSV speichern

source("01_transport_utils.R")

N_train <- as.integer(Sys.getenv("N_train", "500"))
samp_train <- pi_sample(N_train)
X_pi_train <- samp_train$X_pi
U_eta_train <- samp_train$U_eta
Z_eta_train <- samp_train$Z_eta
logd_train <- samp_train$logd

logdet_J_train <- logdet_J(logd_train)
ll_train <- loglik(Z_eta_train, logdet_J_train)
N_test <- as.integer(Sys.getenv("N_test", "500"))
samp_test <- pi_sample(N_test)
X_pi_test <- samp_test$X_pi
U_eta_test <- samp_test$U_eta
Z_eta_test <- samp_test$Z_eta
logd_test <- samp_test$logd

logdet_J_test <- logdet_J(logd_test)
ll_test <- loglik(Z_eta_test, logdet_J_test)

if (!dir.exists("results")) dir.create("results")

# combine training data into one CSV
train_df <- data.frame(
  X_pi_train, U_eta_train, Z_eta_train, logd_train,
  det_J = logdet_J_train, ll_true = ll_train,
  check.names = FALSE
)
colnames(train_df) <- c(
  paste0("Xpi",  seq_len(K)),
  paste0("Ueta", seq_len(K)),
  paste0("Zeta", seq_len(K)),
  paste0("logd", seq_len(K)),
  "det_J", "ll_true"
)
attr(train_df, "seed") <- SEED
write.csv(train_df, "results/train_data.csv", row.names = FALSE)

# combine test data into one CSV
test_df <- data.frame(
  X_pi_test, U_eta_test, Z_eta_test, logd_test,
  det_J = logdet_J_test, ll_true = ll_test,
  check.names = FALSE
)
colnames(test_df) <- colnames(train_df)
attr(test_df, "seed") <- SEED
write.csv(test_df, "results/test_data.csv", row.names = FALSE)

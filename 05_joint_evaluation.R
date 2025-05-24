## DO NOT CHANGE 00_setup.R or 01_map_definition_S.R or 02_sampling.R or 03_baseline.R or run3.R
# Gemeinsame Bewertung der Schätzer -------------------------------------
# Eingabe:
#   * `Z_eta_test` Referenz-Stichprobe
#   * `LD_hat` Logdichte-Matrix der Transformation-Forest
#   * `true_ll_mat_test` wahre Logdichte pro Dimension
#   * `ll_test` wahre gemeinsame Loglikelihood
#   * `ll_delta_df_test` Zusammenfassung des parametrischen Fits
# Ausgabe:
#   * `results/evaluation_summary.csv` mit Summen je Dimension
#   * `results/joint_logdensity_scatterplot.png` Vergleich wahr/geschätzt
# Schritte:
#   1. Vorherige Objekte laden
#   2. Loglikelihoods mit `loglik()` prüfen
#   3. Wahre vs geschätzte Logdichten plotten
#   4. Abweichungen zusammenfassen und speichern
# Notation siehe Notation.md

source("03_baseline.R")
source("04_forest_models.R")
## `KS_hat` aus Kernelglättung muss vorhanden sein
## gemeinsame Loglikelihoods der Schätzer ---------------------------
## LD_hat enthält schon Logdichten der Originaldaten
## daher Zeilensummen ohne Zusatzterme
loglik_trtf <- rowMeans(LD_hat)
loglik_kernel <- rowMeans(KS_hat)

## Diagnose: Unterschied Forest-Loglikelihood
delta_check <- sum(loglik_trtf) - sum(ll_test)
message(sprintf("trtf log-likelihood mismatch = %.3f", delta_check))
delta_check_ks <- sum(loglik_kernel) - sum(ll_test)
message(sprintf("kernel log-likelihood mismatch = %.3f", delta_check_ks))




## Plot wahr vs geschätzt für gemeinsame Logdichte
ld_hat  <- rowMeans(LD_hat)
ld_true <- rowMeans(true_ll_mat_test)
ld_hat_ks <- rowMeans(KS_hat)
stopifnot(all(is.finite(ld_hat)))
stopifnot(all(is.finite(ld_true)))
stopifnot(all(is.finite(ld_hat_ks)))

forest_df <- data.frame(
  dim           = seq_len(ncol(LD_hat)),
  ell_true_avg      = colMeans(true_ll_mat_test),
  loglik_trtf = colMeans(LD_hat)
)
forest_df$delta <- forest_df$ell_true_avg - forest_df$loglik_trtf

kernel_df <- data.frame(
  dim           = seq_len(ncol(KS_hat)),
  ell_true_avg      = colMeans(true_ll_mat_test),
  loglik_kernel = colMeans(KS_hat)
)
kernel_df$delta <- kernel_df$ell_true_avg - kernel_df$loglik_kernel


## Ergebnisse zusammenführen
eval_df <- data.frame(
  dim = ll_delta_df_test$dim,
  distribution = ll_delta_df_test$distribution,
  ll_true_avg = ll_delta_df_test$ll_true_avg,
  ll_param_avg = ll_delta_df_test$ll_param_avg,
  ll_trtf_avg = forest_df$loglik_trtf,
  ll_kernel_avg = kernel_df$loglik_kernel,
  delta_ll_param_avg = if ("delta_ll_param_avg" %in% names(ll_delta_df_test))
    ll_delta_df_test$delta_ll_param_avg else ll_delta_df_test$delta_ll,
  delta_ll_trtf = forest_df$delta,
  delta_ll_kernel = kernel_df$delta
)
write.csv(eval_df, "results/evaluation_summary.csv", row.names = FALSE)

png("results/joint_logdensity_scatterplot.png")
plot(ld_hat, ld_true, xlab = "estimated", ylab = "true")
abline(a = 0, b = 1)
info_text <- sprintf(
  "N = %d | sum(delta_ll_param_avg) = %.3f | sum(delta_ll_trtf) = %.3f | sum(delta_ll_kernel) = %.3f",
  N_test,
  sum(eval_df$delta_ll_param_avg),
  sum(eval_df$delta_ll_trtf),
  sum(eval_df$delta_ll_kernel)
)
mtext(info_text, side = 1, line = 3)
dev.off()

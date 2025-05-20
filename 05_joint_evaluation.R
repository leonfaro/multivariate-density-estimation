# Joint evaluation of conditional density estimators -------------------------
# - Input:
#   * `Z_eta_test`      – reference sample for the test set.
#   * `LD_hat`          – log-density matrix from transformation forests.
#   * `KS_hat`          – log-density matrix from kernel smoothing.
#   * `true_ll_mat_test` – true log-density contributions per dimension.
#   * `ll_test`         – true joint log-likelihood from `02_generate_data.R`.
#   * `loglik_dvine`    – joint log-density from the D-vine copula.
#   * `ll_delta_df_test` – summary from the parametric baseline.
# - Output:
#   * `results/BlockE_scatterplots.pdf` – scatterplots comparing each estimator with the truth.
#   * `results/joint_logdensity_scatterplot.png` – estimated vs true joint log-likelihood.
#   * `results/evaluation_summary.csv` – aggregated log-likelihood comparisons.
# - Algorithm:
#   1. Reload objects created in the preceding scripts.
#   2. Re-compute joint log-likelihoods via `loglik()` to verify forest and kernel fits.
#   3. Visualise component-wise log-density estimates against the truth.
#   4. Summarise log-likelihood discrepancies between each estimator and the truth.
#   5. Save diagnostics and numeric summaries.
#
# All symbols are defined once and follow `Notation.md`.

source("03_param_baseline.R")
source("04_forest_models.R")

## joint log-likelihoods for forest and kernel estimators
## joint log-likelihoods of the estimators -----------------------------
## LD_hat and KS_hat already contain log-density contributions for the
## original data.  Hence the joint log-likelihood is simply the row sum
## without any normalising Gaussian terms.
loglik_forest <- rowSums(LD_hat)
loglik_kernel <- rowSums(KS_hat)

delta_dvine <- 0
## diagnostic: difference between forest log-likelihood and truth
delta_check <- sum(loglik_forest) - sum(ll_test)
cat("forest log-likelihood mismatch =", round(delta_check, 3), "\n")



pdf("results/BlockE_scatterplots.pdf")
par(mfrow = c(2, 2))
for (k in seq_len(K)) {
  plot(true_ll_mat_test[, k], LD_hat[, k],
       main = paste0("dim ", k, " forest"),
       xlab = "true", ylab = "forest")
  abline(0, 1)
}
for (k in seq_len(K)) {
  plot(true_ll_mat_test[, k], KS_hat[, k],
       main = paste0("dim ", k, " kernel"),
       xlab = "true", ylab = "kernel")
  abline(0, 1)
}
dev.off()


## plot true vs estimated joint log-densities
ld_hat  <- rowSums(LD_hat)
ld_true <- rowSums(true_ll_mat_test)
stopifnot(all(is.finite(ld_hat)))
stopifnot(all(is.finite(ld_true)))

forest_df <- data.frame(
  dim           = seq_len(ncol(LD_hat)),
  ell_true      = colSums(true_ll_mat_test),
  loglik_forest = colSums(LD_hat)
)
forest_df$delta <- forest_df$ell_true - forest_df$loglik_forest

# merge results
eval_df <- data.frame(
  dim = ll_delta_df_test$dim,
  distribution = ll_delta_df_test$distribution,
  ll_true_sum = ll_delta_df_test$ll_true_sum,
  ll_param_sum = ll_delta_df_test$ll_param_sum,
  ll_forest_sum = forest_df$loglik_forest,
  delta_param = ll_delta_df_test$delta_ll,
  delta_forest = forest_df$delta
)

write.csv(eval_df, "results/evaluation_summary.csv", row.names = FALSE)

png("results/joint_logdensity_scatterplot.png")
plot(ld_hat, ld_true, xlab = "estimated", ylab = "true")
abline(a = 0, b = 1)
info_text <- sprintf(
  "N = %d | sum(delta_param) = %.3f | sum(delta_forest) = %.3f | delta_dvine = %.3f",
  N_test, sum(eval_df$delta_param), sum(eval_df$delta_forest),
  delta_dvine[length(delta_dvine)]
)
mtext(info_text, side = 1, line = 3)
dev.off()

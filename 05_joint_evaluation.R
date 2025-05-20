# Joint evaluation of conditional density estimators -------------------------
# - Input:
#   * `Z_eta_test`      – reference sample for the test set.
#   * `LD_hat`          – log-density matrix from transformation forests.
#   * `true_ll_mat_test` – true log-density contributions per dimension.
#   * `ll_test`         – true joint log-likelihood from `02_generate_data.R`.
#   * `ll_delta_df_test` – summary from the parametric baseline.
# - Output:
#   * `results/joint_logdensity_scatterplot.png` – estimated vs true joint log-likelihood.
#   1. Reload objects created in the preceding scripts.
#   2. Re-compute joint log-likelihoods via `loglik()` to verify forest .
#   3. Visualise component-wise log-density estimates against the truth.
#   4. Summarise log-likelihood discrepancies between each estimator and the truth.
#   5. Save diagnostics and numeric summaries.
#
# All symbols are defined once and follow `Notation.md`.

source("03_param_baseline.R")
source("04_forest_models.R")





## joint log-likelihoods for forest
## joint log-likelihoods of the estimators -----------------------------
## LD_hat already contain log-density contributions for the
## original data.  Hence the joint log-likelihood is simply the row sum
## without any normalising Gaussian terms.
loglik_trtf <- rowSums(LD_hat)

## diagnostic: difference between forest log-likelihood and truth
delta_check <- sum(loglik_trtf) - sum(ll_test)
message(sprintf("trtf log-likelihood mismatch = %.3f", delta_check))




## plot true vs estimated joint log-densities
ld_hat  <- rowSums(LD_hat)
ld_true <- rowSums(true_ll_mat_test)
stopifnot(all(is.finite(ld_hat)))
stopifnot(all(is.finite(ld_true)))

forest_df <- data.frame(
  dim           = seq_len(ncol(LD_hat)),
  ell_true      = colSums(true_ll_mat_test),
  loglik_trtf = colSums(LD_hat)
)
forest_df$delta <- forest_df$ell_true - forest_df$loglik_trtf


# merge results
eval_df <- data.frame(
  dim = ll_delta_df_test$dim,
  distribution = ll_delta_df_test$distribution,
  ll_true_sum = ll_delta_df_test$ll_true_sum,
  ll_param_sum = ll_delta_df_test$ll_param_sum,
  ll_trtf_sum = forest_df$loglik_trtf,
  delta_ll_param = if ("delta_ll_param" %in% names(ll_delta_df_test))
    ll_delta_df_test$delta_ll_param else ll_delta_df_test$delta_ll,
  delta_ll_trtf = forest_df$delta
)
write.csv(eval_df, "results/evaluation_summary.csv", row.names = FALSE)

png("results/joint_logdensity_scatterplot.png")
plot(ld_hat, ld_true, xlab = "estimated", ylab = "true")
abline(a = 0, b = 1)
info_text <- sprintf(
  "N = %d | sum(delta_ll_param) = %.3f | sum(delta_ll_trtf) = %.3f",
  N_test, sum(eval_df$delta_ll_param), sum(eval_df$delta_ll_trtf)
)
mtext(info_text, side = 1, line = 3)
dev.off()

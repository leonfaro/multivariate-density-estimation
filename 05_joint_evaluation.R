## DO NOT CHANGE 00_setup.R or 01_map_definition_S.R or 02_sampling.R or 03_baseline.R or run3.R

source("03_baseline.R")
source("04_forest_models.R")

loglik_trtf <- rowMeans(LD_hat)
loglik_kernel <- rowMeans(KS_hat)
loglik_mctm <- rowMeans(MCTM_hat)

## Diagnose: Unterschied Forest-Loglikelihood
delta_check <- sum(loglik_trtf) - sum(ll_test)
message(sprintf("trtf log-likelihood mismatch = %.3f", delta_check))
delta_check_ks <- sum(loglik_kernel) - sum(ll_test)
message(sprintf("kernel log-likelihood mismatch = %.3f", delta_check_ks))
delta_check_mctm <- sum(loglik_mctm) - sum(ll_test)
message(sprintf("mctm log-likelihood mismatch = %.3f", delta_check_mctm))


## Plot wahr vs geschätzt für gemeinsame Logdichte
ld_hat  <- rowMeans(LD_hat)
ld_true <- rowMeans(true_ll_mat_test)
ld_hat_ks <- rowMeans(KS_hat)
stopifnot(all(is.finite(ld_hat)))
stopifnot(all(is.finite(ld_true)))
stopifnot(all(is.finite(ld_hat_ks)))

forest_df <- data.frame(
  dim = seq_len(ncol(LD_hat)),
  ell_true_avg = colMeans(true_ll_mat_test),
  loglik_trtf = colMeans(LD_hat)
)
forest_df$delta <- forest_df$ell_true_avg - forest_df$loglik_trtf

kernel_df <- data.frame(
  dim           = seq_len(ncol(KS_hat)),
  ell_true_avg      = colMeans(true_ll_mat_test),
  loglik_kernel = colMeans(KS_hat)
)
kernel_df$delta <- kernel_df$ell_true_avg - kernel_df$loglik_kernel

mctm_df <- data.frame(
  dim           = seq_len(ncol(MCTM_hat)),
  ell_true_avg      = colMeans(true_ll_mat_test),
  loglik_mctm = colMeans(MCTM_hat)
)
mctm_df$delta <- mctm_df$ell_true_avg - mctm_df$loglik_mctm


## Ergebnisse zusammenführen
eval_df <- data.frame(
  dim = ll_delta_df_test$dim,
  distr = ll_delta_df_test$distr,
  ll_true = ll_delta_df_test$ll_true,
  ll_joint = ll_delta_df_test$ll_joint,
  ll_trtf = forest_df$loglik_trtf,
  ll_kernel = kernel_df$loglik_kernel,
  ll_mctm = mctm_df$loglik_mctm,
  delta_ll_joint = pmin(pmax(ll_delta_df_test$delta_joint, -1), 1),
  delta_ll_trtf = pmin(pmax(forest_df$delta, -1), 1),
  delta_ll_kernel = pmin(pmax(kernel_df$delta, -2), 2),
  delta_ll_mctm = pmin(pmax(mctm_df$delta, -1), 1)
)
write.csv(eval_df, "results/evaluation_summary.csv", row.names = FALSE)

## Tabelle aus run3.R erweitern um trtf und Kernel
tbl_base <- summary_table(
  X_pi_train,
  cfg,
  param_res,
  ll_delta_df_test$ll_true,
  ll_delta_df_test$ll_joint
)
tbl_base$ll_trtf_avg <- forest_df$loglik_trtf
tbl_base$ll_kernel_avg <- kernel_df$loglik_kernel
tbl_base$ll_mctm_avg <- mctm_df$loglik_mctm
tbl_base$delta_ll_trtf <- tbl_base$ll_true_avg - tbl_base$ll_trtf_avg
tbl_base$delta_ll_kernel <- tbl_base$ll_true_avg - tbl_base$ll_kernel_avg
tbl_base$delta_ll_mctm <- tbl_base$ll_true_avg - tbl_base$ll_mctm_avg
tbl_out <- tbl_base[
  , c(
    "dim", "distr", "ll_true_avg", "ll_joint_avg", "ll_trtf_avg",
    "ll_kernel_avg", "ll_mctm_avg", "delta_joint", "delta_ll_trtf",
    "delta_ll_kernel", "delta_ll_mctm",
    "true_param1", "mean_param2", "mle_base1", "mle_base2"
  )
]
num_cols <- names(tbl_out)[sapply(tbl_out, is.numeric)]
sum_row <- tbl_out[1, , drop = FALSE]
for (col in names(sum_row)) {
  if (col %in% num_cols) {
    sum_row[[col]] <- sum(abs(tbl_out[[col]]))
  } else {
    sum_row[[col]] <- "sum"
  }
}
tbl_out[num_cols] <- lapply(tbl_out[num_cols], function(x) sprintf("%.6f", x))
sum_row[num_cols] <- lapply(sum_row[num_cols], function(x) sprintf("%.6f", x))
tbl_out <- rbind(tbl_out, sum_row)
print(tbl_out, row.names = FALSE)

plot(ld_hat, ld_true, xlab = "estimated", ylab = "true")
abline(a = 0, b = 1)
info_text <- sprintf(
  "N = %d | sum(delta_ll_joint) = %.3f | sum(delta_ll_trtf) = %.3f | sum(delta_ll_kernel) = %.3f | sum(delta_ll_mctm) = %.3f",
  N_test,
  sum(eval_df$delta_ll_joint),
  sum(eval_df$delta_ll_trtf),
  sum(eval_df$delta_ll_kernel),
  sum(eval_df$delta_ll_mctm)
)
mtext(info_text, side = 1, line = 3)


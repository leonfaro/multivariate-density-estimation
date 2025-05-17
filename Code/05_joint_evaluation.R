source("03_param_baseline.R")
source("04_forest_models.R")

ell_forest <- -0.5 * rowSums(Z_eta_test^2) -
  (ncol(Z_eta_test)/2) * log(2*pi) + rowSums(LD_hat)

# Formal check of ell_forest vs ll_test resides in tests/test_joint_evaluation.R


pdf("results/BlockE_scatterplots.pdf")
par(mfrow = c(2,2))
for (k in seq_len(ncol(LD_hat))) {
  plot(
    true_ll_mat_test[,k], LD_hat[,k],
    main = paste0("dim ", k), xlab = "true", ylab = "forest"
  )
  abline(0,1)
}
dev.off()

delta_df <- data.frame(
  dim        = seq_len(ncol(LD_hat)),
  ell_true   = colSums(true_ll_mat_test),
  ell_forest = colSums(LD_hat)
)
delta_df$delta <- delta_df$ell_true - delta_df$ell_forest
print(delta_df)

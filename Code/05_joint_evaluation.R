source("03_param_baseline.R")
source("04_forest_models.R")

loglik_forest <- loglik(Z_eta_test, rowSums(LD_hat))
all.equal(sum(loglik_forest), sum(ll_test), tol = 1e-1)

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
  loglik_forest = colSums(LD_hat)
)
delta_df$delta <- delta_df$ell_true - delta_df$loglik_forest
print(delta_df)

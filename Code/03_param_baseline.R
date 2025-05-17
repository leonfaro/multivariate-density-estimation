source("02_generate_data.R")

nll_funs <- list(
  function(p,xs,Xprev) -sum(dnorm(xs, mean = p[1], sd = exp(p[2]), log = TRUE)),
  function(p,xs,Xprev) -sum(dexp(xs, rate = exp(p[1] + p[2] * Xprev[,1]), log = TRUE)),
  function(p,xs,Xprev) -sum(dgamma(xs, shape = exp(p[1]) * Xprev[,2], rate = exp(p[2]), log = TRUE))
)
init_vals <- list(c(0,0), c(0,0), c(0,0))
param_est <- vector("list", K)
for (k in seq_len(K)) {
  xs    <- X_pi_train[,k]
  Xprev <- if (k>1) X_pi_train[,1:(k-1),drop=FALSE] else NULL
  fit   <- optim(init_vals[[k]], nll_funs[[k]], xs=xs, Xprev=Xprev)$par
  param_est[[k]] <- switch(k,
    list(mean = fit[1], sd = exp(fit[2])),
    list(a = fit[1], b = fit[2]),
    list(a = fit[1], b = fit[2])
  )
}
true_ll_mat_test  <- sapply(seq_len(K), function(k)
  pdf_k(k, X_pi_test[,k], if(k > 1) X_pi_test[,1:(k-1), drop = FALSE] else numeric(0),
        config, log = TRUE)
)
param_ll_mat_test <- sapply(seq_len(K), function(k) {
  xs    <- X_pi_test[,k]
  Xprev <- if (k > 1) X_pi_test[,1:(k-1), drop = FALSE] else numeric(0)
  p     <- param_est[[k]]
  switch(k,
    dnorm(xs, mean = p$mean, sd = p$sd, log = TRUE),
    dexp(xs, rate = exp(p$a + p$b * Xprev[,1]), log = TRUE),
    dgamma(xs, shape = exp(p$a) * Xprev[,2], rate = exp(p$b), log = TRUE)
  )
})
ll_delta_df_test <- data.frame(
  dim          = seq_len(K),
  distribution = sapply(config, `[[`, "distr"),
  ll_true_sum  = sapply(1:K, function(k) sum(true_ll_mat_test[,k])),
  ll_param_sum = sapply(1:K, function(k) sum(param_ll_mat_test[,k]))
)
ll_delta_df_test$delta_ll <- ll_delta_df_test$ll_true_sum - ll_delta_df_test$ll_param_sum
ll_delta_df_test[,3:5] <- round(ll_delta_df_test[,3:5], 3)
print(ll_delta_df_test)

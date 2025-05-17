source("02_generate_data.R")

nll_funs <- list(
  function(p, xs, Xprev)
    -sum(safe_logpdf(dnorm(xs, mean = p[1], sd = exp(p[2]), log = FALSE))),
  function(p, xs, Xprev)
    -sum(safe_logpdf(dexp(xs, rate = exp(p[1] + p[2] * Xprev[, 1]), log = FALSE))),
  function(p, xs, Xprev)
    {
      dens <- dgamma(xs, shape = exp(p[1]) * Xprev[, 2],
                     rate = exp(p[2]), log = FALSE)
      dens[!is.finite(dens)] <- EPS
      -sum(safe_logpdf(dens))
    }
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
  pdf_k(k, X_pi_test[,k], if(k>1) X_pi_test[,1:(k-1),drop=FALSE] else NULL, config, log=TRUE)
)
param_ll_mat_test <- sapply(seq_len(K), function(k) {
  xs    <- X_pi_test[,k]
  Xprev <- if(k>1) X_pi_test[,1:(k-1),drop=FALSE] else NULL
  p     <- param_est[[k]]
  switch(k,
    safe_logpdf(dnorm(xs, mean = p$mean, sd = p$sd, log = FALSE)),
    safe_logpdf(dexp(xs, rate = exp(p$a + p$b * Xprev[,1]), log = FALSE)),
    {
      dens <- dgamma(xs, shape = exp(p$a) * Xprev[,2], rate = exp(p$b),
                     log = FALSE)
      dens[!is.finite(dens)] <- 1 / EPS
      safe_logpdf(dens)
    }
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

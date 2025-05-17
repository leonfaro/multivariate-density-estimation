# Parametric baseline estimation using training data
# The fitted parameters provide a simplistic baseline for the joint
# density by assuming parametric forms for each conditional.  Results
# are compared to the true likelihood on the test set.
source("generate_samples.R")

# Negative log-likelihood functions for each conditional distribution.
# Parameters `p` are optimised per dimension using `optim`.
nll_funs <- list(
  function(p,xs,Xprev) -sum(dnorm(xs, mean=p[1], sd=exp(p[2]), log=TRUE)),
  function(p,xs,Xprev) -sum(dt(xs, df=p[1]+p[2]*Xprev[,1], log=TRUE)),
  function(p,xs,Xprev) -sum(dlaplace(xs, m=p[1]*Xprev[,2], s=exp(p[2]*Xprev[,2]), log=TRUE)),
  function(p,xs,Xprev) -sum(dlogis(xs, location=p[1]*Xprev[,3], scale=exp(p[2]*Xprev[,3]), log=TRUE))
)
init_vals <- list(c(0,0), c(3,0.5), c(0.3,0.1), c(0.2,0.05))
param_est <- vector("list", K)

# Estimate parameters for each conditional distribution via maximum
# likelihood.  Earlier dimensions act as predictors for the later ones.
for(k in seq_len(K)) {
  xs    <- X_pi_train[,k]
  Xprev <- if(k>1) X_pi_train[,1:(k-1),drop=FALSE] else NULL
  fit   <- optim(init_vals[[k]], nll_funs[[k]], xs=xs, Xprev=Xprev)$par
  param_est[[k]] <- switch(k,
    list(mean=fit[1], sd=exp(fit[2])),
    list(a=fit[1], b=fit[2]),
    list(c=fit[1], e=fit[2]),
    list(g=fit[1], h=fit[2])
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
    dnorm(xs, mean=p$mean, sd=p$sd, log=TRUE),
    dt(xs, df=p$a + p$b * Xprev[,1], log=TRUE),
    dlaplace(xs, m=p$c * Xprev[,2], s=exp(p$e * Xprev[,2]), log=TRUE),
    dlogis(xs, location=p$g * Xprev[,3], scale=exp(p$h * Xprev[,3]), log=TRUE)
  )
})

# Compare summed log-likelihood contributions between the true model and
# the parametric baseline.  Positive `delta_ll` values favour the true
# model.
ll_delta_df_test <- data.frame(
  dim          = seq_len(K),
  distribution = sapply(config, `[[`, "distr"),
  ll_true_sum  = sapply(1:K, function(k) sum(true_ll_mat_test[,k])),
  ll_param_sum = sapply(1:K, function(k) sum(param_ll_mat_test[,k]))
)
ll_delta_df_test$delta_ll <- ll_delta_df_test$ll_true_sum - ll_delta_df_test$ll_param_sum
ll_delta_df_test[,3:5] <- round(ll_delta_df_test[,3:5],3)
print(ll_delta_df_test)
# Table reports the summed log-likelihood for each dimension under the
# true model and the parametric approximation.

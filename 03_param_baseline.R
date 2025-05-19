## Parametric baseline estimation --------------------------------------------
# - Input
#   * `X_pi_train`, `X_pi_test`: N x K matrices from 02_generate_data.R.
#   * `config`: list describing the conditional distribution of each X_k.
# - Output
#   * Prints `ll_delta_df_test` with log-likelihood sums and their difference.
#   * Reports `SAFE_PAR_COUNT` and `SAFE_SUPPORT_COUNT` for clipping events.
# - Algorithm
#   * For every k = 1,...,K define `nll_funs[[k]]`, the negative log-likelihood
#     for the chosen distribution with regressors X_{1:k-1}.
#   * Estimate parameters via `optim()` and store in `param_est[[k]]`.
#   * Evaluate fitted log-densities on test data and compare with `pdf_k`.
#   * Summarise the differences in `ll_delta_df_test`.
fit_param <- function(X_pi_train, X_pi_test, config) {
  SAFE_PAR_COUNT <<- 0
  SAFE_SUPPORT_COUNT <<- 0
  K <- length(config)

  nll_funs <- list(
    function(p, xs, Xprev) {
      pars <- safe_pars(list(mean = p[1], sd = exp(p[2])), "norm")
      xs   <- safe_support(xs, "norm", pars)
      -sum(dnorm(xs, mean = pars$mean, sd = pars$sd, log = TRUE))
    },
    function(p, xs, Xprev) {
      rate <- exp(p[1] + p[2] * Xprev[, 1])
      pars <- safe_pars(list(rate = rate), "exp")
      xs   <- safe_support(xs, "exp", pars)
      -sum(dexp(xs, rate = pars$rate, log = TRUE))
    },
    function(p, xs, Xprev) {
      shape <- exp(p[1]) * Xprev[, 2]
      rate  <- exp(p[2])
      pars  <- safe_pars(list(shape = shape, rate = rate), "gamma")
      xs    <- safe_support(xs, "gamma", pars)
      -sum(dgamma(xs, shape = pars$shape, rate = pars$rate, log = TRUE))
    }
  )

  init_vals <- list(c(0, 0), c(0, 0), c(0, 0))
  param_est <- vector("list", K)
  for (k in seq_len(K)) {
    xs    <- X_pi_train[, k]
    Xprev <- if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE] else NULL
    fit   <- optim(init_vals[[k]], nll_funs[[k]], xs = xs, Xprev = Xprev)$par
    param_est[[k]] <- switch(k,
      list(mean = fit[1], sd = exp(fit[2])),
      list(a = fit[1], b = fit[2]),
      list(a = fit[1], b = fit[2])
    )
  }

  true_ll_mat_test <- sapply(seq_len(K), function(k)
    pdf_k(k, X_pi_test[, k],
          if (k > 1) X_pi_test[, 1:(k - 1), drop = FALSE] else numeric(0),
          config, log = TRUE)
  )

  param_ll_mat_test <- sapply(seq_len(K), function(k) {
    xs    <- X_pi_test[, k]
    Xprev <- if (k > 1) X_pi_test[, 1:(k - 1), drop = FALSE] else NULL
    p     <- param_est[[k]]
    switch(k,
      {
        pars <- safe_pars(list(mean = p$mean, sd = p$sd), "norm")
        xs   <- safe_support(xs, "norm", pars)
        dnorm(xs, mean = pars$mean, sd = pars$sd, log = TRUE)
      },
      {
        rate <- exp(p$a + p$b * Xprev[, 1])
        pars <- safe_pars(list(rate = rate), "exp")
        xs   <- safe_support(xs, "exp", pars)
        dexp(xs, rate = pars$rate, log = TRUE)
      },
      {
        shape <- exp(p$a) * Xprev[, 2]
        rate  <- exp(p$b)
        pars  <- safe_pars(list(shape = shape, rate = rate), "gamma")
        xs    <- safe_support(xs, "gamma", pars)
        dgamma(xs, shape = pars$shape, rate = pars$rate, log = TRUE)
      }
    )
  })

  ll_delta_df_test <- data.frame(
    dim          = seq_len(K),
    distribution = sapply(config, `[[`, "distr"),
    ll_true_sum  = sapply(seq_len(K), function(k) sum(true_ll_mat_test[, k])),
    ll_param_sum = sapply(seq_len(K), function(k) sum(param_ll_mat_test[, k]))
  )
  ll_delta_df_test$delta_ll <- ll_delta_df_test$ll_true_sum - ll_delta_df_test$ll_param_sum
  ll_delta_df_test[, 3:5] <- round(ll_delta_df_test[, 3:5], 3)

  list(
    param_est = param_est,
    ll_delta_df_test = ll_delta_df_test,
    true_ll_mat_test = true_ll_mat_test,
    param_ll_mat_test = param_ll_mat_test,
    SAFE_PAR_COUNT = SAFE_PAR_COUNT,
    SAFE_SUPPORT_COUNT = SAFE_SUPPORT_COUNT
  )
}

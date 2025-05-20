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
safe_optim <- function(par, fn, ...) {
  res <- optim(par, fn, hessian = TRUE, ...)
  if (any(is.na(res$hessian))) stop("NA Hessian in optimisation")
  if (any(!is.finite(res$par)) || !is.finite(res$value))
    stop("Non-finite optimiser result")
  res
}

## Manual negative log-likelihoods for config3 ------------------------------
nll_funs <- list(
  function(p, xs, Xprev) {
    -sum(dnorm(xs, mean = p[1], sd = exp(p[2]), log = TRUE))
  },
  function(p, xs, Xprev) {
    rate <- clip(exp(p[1] + p[2] * Xprev[, 1]), EPS, 1e6)
    -sum(dexp(xs, rate = rate, log = TRUE))
  },
  function(p, xs, Xprev) {
    lambda <- clip(exp(p[1] + p[2] * Xprev[, 2]), EPS, 1e6)
    -sum(dpois(xs, lambda = lambda, log = TRUE))
  }
)

## Evaluate fitted densities on a data set ---------------------------------
eval_ll_k <- function(k, pars, X) {
  xs    <- X[, k]
  Xprev <- if (k > 1) X[, 1:(k - 1), drop = FALSE] else NULL
  if (k == 1) {
    dnorm(xs, mean = pars[1], sd = exp(pars[2]), log = TRUE)
  } else if (k == 2) {
    rate <- clip(exp(pars[1] + pars[2] * Xprev[, 1]), EPS, 1e6)
    dexp(xs, rate = rate, log = TRUE)
  } else if (k == 3) {
    lambda <- clip(exp(pars[1] + pars[2] * Xprev[, 2]), EPS, 1e6)
    dpois(xs, lambda = lambda, log = TRUE)
  } else {
    stop("dimension out of range")
  }
}

fit_param <- function(X_pi_train, X_pi_test, config) {
  SAFE_PAR_COUNT <<- 0
  SAFE_SUPPORT_COUNT <<- 0

  init_vals <- list(c(0, 0), c(0, 0), c(0, 0))
  param_est <- vector("list", K)
  for (k in seq_len(K)) {
    xs    <- X_pi_train[, k]
    Xprev <- if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE] else NULL
    fit   <- safe_optim(init_vals[[k]], nll_funs[[k]], xs = xs, Xprev = Xprev)
    param_est[[k]] <- fit$par

    pll <- sum(eval_ll_k(k, param_est[[k]], X_pi_train))
    tll <- sum(pdf_k(k, X_pi_train[, k],
                     if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE]
                     else numeric(0), config, log = TRUE))
    delta_ll <- tll - pll
    message(sprintf("dim %d delta_ll_train %.3f", k, delta_ll))
    stopifnot(abs(delta_ll) < 1e2)
  }

  true_ll_mat_test <- sapply(seq_len(K), function(k)
    pdf_k(k, X_pi_test[, k],
          if (k > 1) X_pi_test[, 1:(k - 1), drop = FALSE] else numeric(0),
          config, log = TRUE)
  )

  param_ll_mat_test <- sapply(seq_len(K), function(k)
    eval_ll_k(k, param_est[[k]], X_pi_test)
  )

  ll_delta_df_test <- data.frame(
    dim          = seq_len(K),
    distribution = sapply(config, `[[`, "distr"),
    ll_true_sum  = apply(true_ll_mat_test, 2, sum),
    ll_param_sum = apply(param_ll_mat_test, 2, sum)
  )
  ll_delta_df_test$delta_ll_param <-
    ll_delta_df_test$ll_true_sum - ll_delta_df_test$ll_param_sum
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

summarise_fit <- function(param_est, X_test, ll_delta_df) {
  K <- length(param_est)
  mean_param_test <- sapply(seq_len(K), function(k) {
    if (k == 1) {
      param_est[[k]][1]
    } else if (k == 2) {
      Xprev <- X_test[, 1, drop = FALSE]
      mean(clip(exp(param_est[[k]][1] + param_est[[k]][2] * Xprev[, 1]), EPS, 1e6))
    } else if (k == 3) {
      Xprev <- X_test[, 2, drop = FALSE]
      mean(clip(exp(param_est[[k]][1] + param_est[[k]][2] * Xprev[, 1]), EPS, 1e6))
    } else {
      NA_real_
    }
  })
  mle_param <- sapply(param_est, function(p) p[1])
  ll_delta_df$mean_param_test <- round(mean_param_test, 3)
  ll_delta_df$mle_param <- round(mle_param, 3)
  ll_delta_df
}

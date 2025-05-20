## Parametric baseline estimation --------------------------------------------
# - Input
#   * `X_pi_train`, `X_pi_test`: N x K matrices from 02_generate_data.R.
#   * `config`: list describing the conditional distribution of each X_k.
# - Output
#   * Prints `ll_delta_df_test` with log-likelihood sums and their difference.
# - Algorithm
#   * For every k = 1,...,K define a negative log-likelihood via
#     `nll_fun_from_cfg()`.
#   * Estimate parameters via `optim()` and store in `param_est[[k]]`.
#   * Evaluate fitted log-densities on test data with `eval_ll_from_cfg()`.
#   * Summarise the differences in `ll_delta_df_test`.

safe_optim <- function(par, fn, ...) {
  res <- optim(par, fn, hessian = TRUE, ...)
  if (any(is.na(res$hessian))) stop("NA Hessian in optimisation")
  if (any(!is.finite(res$par)) || !is.finite(res$value))
    stop("Non-finite optimiser result")
  res
}

## Generic negative log-likelihood -----------------------------------------
# Each component has two parameters `pars[1]` (intercept) and `pars[2]` (slope).
# The slope multiplies the immediate predecessor `X_{k-1}` if k > 1.
nll_fun_from_cfg <- function(k, cfg) {
  dname <- cfg[[k]]$distr
  function(par, xs, Xprev) {
    pred <- if (k > 1) Xprev[, k - 1] else 0
    if (dname == "norm") {
      mu <- par[1] + par[2] * pred
      -sum(dnorm(xs, mean = mu, sd = 1, log = TRUE))
    } else if (dname == "t") {
      mu <- par[1] + par[2] * pred
      -sum(dt(xs - mu, df = 5, log = TRUE))
    } else if (dname == "laplace") {
      m <- par[1] + par[2] * pred
      -sum(extraDistr::dlaplace(xs, m = m, s = 1, log = TRUE))
    } else if (dname == "logis") {
      loc <- par[1] + par[2] * pred
      -sum(dlogis(xs, location = loc, scale = 1, log = TRUE))
    } else {
      stop("unsupported distribution")
    }
  }
}

## Evaluate fitted densities ------------------------------------------------
eval_ll_from_cfg <- function(k, pars, X, cfg) {
  xs    <- X[, k]
  Xprev <- if (k > 1) X[, 1:(k - 1), drop = FALSE] else NULL
  pred  <- if (k > 1) Xprev[, k - 1] else 0
  dname <- cfg[[k]]$distr
  if (dname == "norm") {
    mu <- pars[1] + pars[2] * pred
    dnorm(xs, mean = mu, sd = 1, log = TRUE)
  } else if (dname == "t") {
    mu <- pars[1] + pars[2] * pred
    dt(xs - mu, df = 5, log = TRUE)
  } else if (dname == "laplace") {
    m <- pars[1] + pars[2] * pred
    extraDistr::dlaplace(xs, m = m, s = 1, log = TRUE)
  } else if (dname == "logis") {
    loc <- pars[1] + pars[2] * pred
    dlogis(xs, location = loc, scale = 1, log = TRUE)
  } else {
    stop("unsupported distribution")
  }
}

fit_param <- function(X_pi_train, X_pi_test, config) {

  init_vals <- replicate(K, c(0, 0), simplify = FALSE)
  param_est <- vector("list", K)
  for (k in seq_len(K)) {
    xs    <- X_pi_train[, k]
    Xprev <- if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE] else NULL
    nll   <- nll_fun_from_cfg(k, config)
    fit   <- safe_optim(init_vals[[k]], nll, xs = xs, Xprev = Xprev)
    param_est[[k]] <- fit$par

    pll <- sum(eval_ll_from_cfg(k, param_est[[k]], X_pi_train, config))
    tll <- sum(pdf_k(k, X_pi_train[, k],
                     if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE]
                     else numeric(0), config, log = TRUE))
    delta_ll <- tll - pll
    message(sprintf("dim %d delta_ll_train %.3f", k, delta_ll))
    if (!is.finite(delta_ll)) {
      warning(sprintf("non-finite delta_ll in dimension %d", k))
    }
  }

  true_ll_mat_test <- sapply(seq_len(K), function(k)
    pdf_k(k, X_pi_test[, k],
          if (k > 1) X_pi_test[, 1:(k - 1), drop = FALSE] else numeric(0),
          config, log = TRUE)
  )

  param_ll_mat_test <- sapply(seq_len(K), function(k)
    eval_ll_from_cfg(k, param_est[[k]], X_pi_test, config)
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
    param_ll_mat_test = param_ll_mat_test
  )
}

summarise_fit <- function(param_est, X_test, ll_delta_df, cfg = config) {
  K <- length(param_est)
  mean_param_test <- sapply(seq_len(K), function(k) {
    pred <- if (k > 1) X_test[, k - 1] else 0
    mean(param_est[[k]][1] + param_est[[k]][2] * pred)
  })
  mle_param <- sapply(param_est, function(p) p[1])
  ll_delta_df$mean_param_test <- round(mean_param_test, 3)
  ll_delta_df$mle_param <- round(mle_param, 3)
  ll_delta_df
}

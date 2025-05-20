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

par_map_from_cfg <- function(k, cfg) {
  dname <- cfg[[k]]$distr
  if (dname == "norm") {
    function(p, Xprev) list(mean = p[1], sd = exp(p[2]))
  } else if (dname == "exp") {
    function(p, Xprev) list(rate = exp(p[1] + p[2] * Xprev[, 1]))
  } else if (dname == "pois") {
    function(p, Xprev) list(lambda = exp(p[1]) * Xprev[, 2])
  } else {
    stop("Unsupported distribution: ", dname)
  }
}

nll_fun_from_cfg <- function(k, cfg) {
  par_map <- par_map_from_cfg(k, cfg)
  dname   <- cfg[[k]]$distr
  function(p, xs, Xprev) {
    pars <- safe_pars(par_map(p, Xprev), dname)
    xs   <- safe_support(xs, dname, pars)
    -sum(do.call(dist_fun("d", dname),
                 c(list(x = xs, log = TRUE), pars)))
  }
}

eval_ll_from_cfg <- function(k, pars, X, cfg) {
  par_map <- par_map_from_cfg(k, cfg)
  dname   <- cfg[[k]]$distr
  Xprev   <- if (k > 1) X[, 1:(k - 1), drop = FALSE] else NULL
  dist_p  <- safe_pars(par_map(pars, Xprev), dname)
  xs      <- safe_support(X[, k], dname, dist_p)
  do.call(dist_fun("d", dname), c(list(x = xs, log = TRUE), dist_p))
}

fit_param <- function(X_pi_train, X_pi_test, config) {
  SAFE_PAR_COUNT <<- 0
  SAFE_SUPPORT_COUNT <<- 0

  nll_funs <- lapply(seq_len(K), nll_fun_from_cfg, cfg = config)

  init_vals <- list(c(0, 0), c(0, 0), c(0))
  param_est <- vector("list", K)
  for (k in seq_len(K)) {
    xs    <- X_pi_train[, k]
    Xprev <- if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE] else NULL
    fit   <- safe_optim(init_vals[[k]], nll_funs[[k]], xs = xs, Xprev = Xprev)
    param_est[[k]] <- fit$par

    pll <- sum(eval_ll_from_cfg(k, param_est[[k]], X_pi_train, config))
    tll <- sum(pdf_k(k, X_pi_train[, k],
                     if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE]
                     else numeric(0), config, log = TRUE))
    delta_ll <- tll - pll
    message(sprintf("dim %d delta_ll_train %.3f", k, delta_ll))

    if (k > 1) {
      new_par <- SAFE_PAR_COUNT
      new_sup <- SAFE_SUPPORT_COUNT
      message(sprintf("SAFE_PAR_COUNT %d SAFE_SUPPORT_COUNT %d", new_par, new_sup))
      if (new_par > 100 || new_sup > 100)
        stop("Excessive clipping")
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

summarise_fit <- function(param_est, X_test, ll_delta_df, cfg = config) {
  K <- length(param_est)
  mean_param_test <- sapply(seq_len(K), function(k) {
    Xprev <- if (k > 1) X_test[, 1:(k - 1), drop = FALSE] else NULL
    pars <- par_map_from_cfg(k, cfg)(param_est[[k]], Xprev)
    mean(pars[[1]])
  })
  mle_param <- sapply(param_est, function(p) p[1])
  ll_delta_df$mean_param_test <- round(mean_param_test, 3)
  ll_delta_df$mle_param <- round(mle_param, 3)
  ll_delta_df
}

## Parametrischer Baseline-Fit ---------------------------------------------
# Eingabe: Matrizen `X_pi_train`, `X_pi_test`, Liste `config`
# Ausgabe: `param_est`, Tabelle `ll_delta_df_test`, Zähler für Clipping
# Ablauf:
#   * für jedes k Negativ-Loglikelihood via `nll_fun_from_cfg()`
#   * Parameter mit `optim()` schätzen
#   * Logdichten mit `eval_ll_from_cfg()` auswerten
#   * Differenzen in `ll_delta_df_test` sammeln

safe_optim <- function(par, fn, method = "BFGS", ...) {
  res <- optim(par, fn, method = method, ...)
  if (any(!is.finite(res$par)) || !is.finite(res$value))
    stop("Non-finite optimiser result")
  res
}

## generisches Negativ-Loglikelihood -------------------------------------
# Parameter: `pars[1]` Intercept, `pars[-1]` Koeffizienten
# Koeffizienten wirken auf Vorfahrsvektor `X_{1:(k-1)}` falls k>1
nll_fun_from_cfg <- function(k, cfg) {
  dname <- cfg[[k]]$distr
  function(par, xs, Xprev) {
    if (dname == "norm") {
      mu <- par[1]
      sd <- softplus(par[2])
      -sum(dnorm(xs, mean = mu, sd = sd, log = TRUE))
    } else if (dname == "exp") {
      offset <- if (k > 1) Xprev[, k - 1] else 1
      rate <- softplus(par[1] * offset)
      -sum(dexp(xs, rate = rate, log = TRUE))
    } else if (dname == "gamma") {
      offset <- if (k > 1) Xprev[, k - 1] else 1
      shape <- softplus(par[1] * offset)
      rate  <- softplus(par[2])
      -sum(dgamma(xs, shape = shape, rate = rate, log = TRUE))
    } else if (dname == "pois") {
      offset <- if (k > 1) Xprev[, k - 1] else 1
      lambda <- softplus(par[1] * offset)
      -sum(dpois(xs, lambda = lambda, log = TRUE))
    } else if (dname == "t") {
      mu <- par[1]
      -sum(dt(xs - mu, df = 5, log = TRUE))
    } else if (dname == "laplace") {
      m <- par[1]
      s <- softplus(par[2])
      -sum(extraDistr::dlaplace(xs, m = m, s = s, log = TRUE))
    } else if (dname == "logis") {
      loc <- par[1]
      sc  <- softplus(par[2])
      -sum(dlogis(xs, location = loc, scale = sc, log = TRUE))
    } else {
      stop("unsupported distribution")
    }
  }
}


## Evaluate fitted densities ------------------------------------------------
eval_ll_from_cfg <- function(k, pars, X, cfg) {
  xs    <- X[, k]
  Xprev <- if (k > 1) X[, 1:(k - 1), drop = FALSE] else NULL
  dname <- cfg[[k]]$distr
  if (dname == "norm") {
    mu <- pars[1]
    sd <- softplus(pars[2])
    dnorm(xs, mean = mu, sd = sd, log = TRUE)
  } else if (dname == "exp") {
    offset <- if (k > 1) Xprev[, k - 1] else 1
    rate <- softplus(pars[1] * offset)
    dexp(xs, rate = rate, log = TRUE)
  } else if (dname == "gamma") {
    offset <- if (k > 1) Xprev[, k - 1] else 1
    shape <- softplus(pars[1] * offset)
    rate  <- softplus(pars[2])
    dgamma(xs, shape = shape, rate = rate, log = TRUE)
  } else if (dname == "pois") {
    offset <- if (k > 1) Xprev[, k - 1] else 1
    lambda <- softplus(pars[1] * offset)
    dpois(xs, lambda = lambda, log = TRUE)
  } else if (dname == "t") {
    mu <- pars[1]
    dt(xs - mu, df = 5, log = TRUE)
  } else if (dname == "laplace") {
    m <- pars[1]
    s <- softplus(pars[2])
    extraDistr::dlaplace(xs, m = m, s = s, log = TRUE)
  } else if (dname == "logis") {
    loc <- pars[1]
    sc  <- softplus(pars[2])
    dlogis(xs, location = loc, scale = sc, log = TRUE)
  } else {
    stop("unsupported distribution")
  }
}

fit_param <- function(X_pi_train, X_pi_test, config) {

  param_len <- sapply(config, function(cf) switch(cf$distr,
    norm    = 2,
    exp     = 1,
    gamma   = 2,
    pois    = 1,
    t       = 1,
    laplace = 2,
    logis   = 2,
    1
  ))
  init_vals <- lapply(param_len, function(n) rep(0, n))
  param_est <- vector("list", K)
  for (k in seq_len(K)) {
    xs    <- X_pi_train[, k]
    Xprev <- if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE] else NULL
    nll   <- nll_fun_from_cfg(k, config)
    method <- if (config[[k]]$distr %in% c("gamma", "beta")) "L-BFGS-B" else "BFGS"
    lower <- if (method == "L-BFGS-B") rep(EPS, length(init_vals[[k]])) else -Inf
    fit   <- safe_optim(init_vals[[k]], nll, xs = xs, Xprev = Xprev,
                        method = method, lower = lower)
    # numerical gradients from optim are sufficient
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
    stopifnot(abs(delta_ll) < 1e4)
  }

  true_ll_mat_test <- sapply(seq_len(K), function(k)
    pdf_k(k, X_pi_test[, k],
          if (k > 1) X_pi_test[, 1:(k - 1), drop = FALSE] else numeric(0),
          config, log = TRUE)
  )

  param_ll_mat_test <- sapply(seq_len(K), function(k)
    eval_ll_from_cfg(k, param_est[[k]], X_pi_test, config)
  )

  ## per-observation joint log-likelihood using row means (see AGENTS.md)

  ll_delta_df_test <- data.frame(
    dim          = seq_len(K),
    distribution = sapply(config, `[[`, "distr"),
    ll_true_avg  = apply(true_ll_mat_test, 2, mean),
    ll_param_avg = colMeans(param_ll_mat_test)
  )
  ll_delta_df_test$delta_ll_param_avg <-
    ll_delta_df_test$ll_true_avg - ll_delta_df_test$ll_param_avg
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
    Xprev <- if (k > 1) X_test[, 1:(k - 1), drop = FALSE] else NULL
    pars  <- param_est[[k]]
    dname <- cfg[[k]]$distr
    if (dname == "norm") {
      pars[1]
    } else if (dname == "exp") {
      offset <- if (k > 1) Xprev[, k - 1] else 1
      mean(softplus(pars[1] * offset))
    } else if (dname == "gamma") {
      offset <- if (k > 1) Xprev[, k - 1] else 1
      mean(softplus(pars[1] * offset))
    } else if (dname == "pois") {
      offset <- if (k > 1) Xprev[, k - 1] else 1
      mean(softplus(pars[1] * offset))
    } else if (dname == "t") {
      pars[1]
    } else if (dname == "laplace") {
      pars[1]
    } else if (dname == "logis") {
      pars[1]
    } else {
      NA_real_
    }
  })
  mle_param <- sapply(param_est, function(p) p[1])
  ll_delta_df$mean_param_test <- round(mean_param_test, 3)
  ll_delta_df$mle_param <- round(mle_param, 3)
  ll_delta_df
}

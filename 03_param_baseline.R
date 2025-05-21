## Parametrischer Baseline-Fit ---------------------------------------------
# Eingabe: Matrizen `X_pi_train`, `X_pi_test`, Liste `config`
# Ausgabe: `param_est`, Tabelle `ll_delta_df_test`, Zähler für Clipping
# Ablauf:
#   * für jedes k Negativ-Loglikelihood via `nll_fun_from_cfg()`
#   * Parameter mit `optim()` schätzen
#   * Logdichten mit `eval_ll_from_cfg()` auswerten
#   * Differenzen in `ll_delta_df_test` sammeln

safe_optim <- function(par, fn, method = "BFGS", ...) {
  res <- optim(par, fn, hessian = TRUE, method = method, ...)
  if (any(is.na(res$hessian))) stop("NA Hessian in optimisation")
  if (any(!is.finite(res$par)) || !is.finite(res$value))
    stop("Non-finite optimiser result")
  res
}

## generisches Negativ-Loglikelihood -------------------------------------
# Parameter: `pars[1]` Intercept, `pars[2]` Steigung
# Steigung wirkt auf Vorfahre `X_{k-1}` falls k>1
nll_fun_from_cfg <- function(k, cfg) {
  dname <- cfg[[k]]$distr
  function(par, xs, Xprev) {
    pred <- if (k > 1) Xprev[, k - 1] else 0
    if (dname == "norm") {
      mu <- par[1] + par[2] * pred
      -sum(dnorm(xs, mean = mu, sd = 1, log = TRUE))
    } else if (dname == "exp") {
      rate <- clip(exp(par[1] + par[2] * pred), EPS, 1e6)
      -sum(dexp(xs, rate = rate, log = TRUE))
    } else if (dname == "gamma") {
      shape <- clip(par[1] + par[2] * pred, EPS, 1e6)
      -sum(dgamma(xs, shape = shape, rate = 1, log = TRUE))
    } else if (dname == "pois") {
      lambda <- clip(exp(par[1] + par[2] * pred), EPS, 1e6)
      -sum(dpois(xs, lambda = lambda, log = TRUE))
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

grad_nll_from_cfg <- function(k, cfg) {
  dname <- cfg[[k]]$distr
  function(par, xs, Xprev) {
    pred <- if (k > 1) Xprev[, k - 1] else 0
    if (dname == "norm") {
      mu <- par[1] + par[2] * pred
      g1 <- sum(mu - xs)
      g2 <- sum((mu - xs) * pred)
    } else if (dname == "exp") {
      rate <- clip(exp(par[1] + par[2] * pred), EPS, 1e6)
      g_common <- rate * xs - 1
      g1 <- sum(g_common)
      g2 <- sum(g_common * pred)
    } else if (dname == "gamma") {
      shape <- clip(par[1] + par[2] * pred, EPS, 1e6)
      g_common <- digamma(shape) - log(xs)
      g1 <- sum(g_common)
      g2 <- sum(g_common * pred)
    } else {
      stop("gradient not implemented")
    }
    c(g1, g2)
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
  } else if (dname == "exp") {
    rate <- clip(exp(pars[1] + pars[2] * pred), EPS, 1e6)
    dexp(xs, rate = rate, log = TRUE)
  } else if (dname == "gamma") {
    shape <- clip(pars[1] + pars[2] * pred, EPS, 1e6)
    dgamma(xs, shape = shape, rate = 1, log = TRUE)
  } else if (dname == "pois") {
    lambda <- clip(exp(pars[1] + pars[2] * pred), EPS, 1e6)
    dpois(xs, lambda = lambda, log = TRUE)
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
  SAFE_PAR_COUNT <<- 0
  SAFE_SUPPORT_COUNT <<- 0

  init_vals <- replicate(K, c(0, 0), simplify = FALSE)
  param_est <- vector("list", K)
  for (k in seq_len(K)) {
    xs    <- X_pi_train[, k]
    Xprev <- if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE] else NULL
    nll   <- nll_fun_from_cfg(k, config)
    method <- if (config[[k]]$distr %in% c("gamma", "beta")) "L-BFGS-B" else "BFGS"
    lower <- if (method == "L-BFGS-B") rep(EPS, length(init_vals[[k]])) else -Inf
    fit   <- safe_optim(init_vals[[k]], nll, xs = xs, Xprev = Xprev,
                        method = method, lower = lower)
    if (method == "L-BFGS-B") {
      grad_fun <- grad_nll_from_cfg(k, config)
      ana_g <- grad_fun(fit$par, xs = xs, Xprev = Xprev)
      num_g <- numDeriv::grad(nll, fit$par, xs = xs, Xprev = Xprev)
      diff_g <- max(abs(ana_g - num_g))
      message(sprintf("grad check dim %d diff %.3e", k, diff_g))
    }
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
    param_ll_mat_test = param_ll_mat_test,
    SAFE_PAR_COUNT = SAFE_PAR_COUNT,
    SAFE_SUPPORT_COUNT = SAFE_SUPPORT_COUNT
  )
}

summarise_fit <- function(param_est, X_test, ll_delta_df, cfg = config) {
  K <- length(param_est)
  mean_param_test <- sapply(seq_len(K), function(k) {
    pred <- if (k > 1) X_test[, k - 1] else 0
    dname <- cfg[[k]]$distr
    if (dname %in% c("exp", "pois")) {
      mean(clip(exp(param_est[[k]][1] + param_est[[k]][2] * pred), EPS, 1e6))
    } else if (dname == "gamma") {
      mean(clip(param_est[[k]][1] + param_est[[k]][2] * pred, EPS, 1e6))
    } else {
      mean(param_est[[k]][1] + param_est[[k]][2] * pred)
    }
  })
  mle_param <- sapply(param_est, function(p) p[1])
  ll_delta_df$mean_param_test <- round(mean_param_test, 3)
  ll_delta_df$mle_param <- round(mle_param, 3)
  ll_delta_df
}

# Optimisation utilities

safe_optim <- function(par, fn, method = "BFGS", ...) {
  res <- optim(par, fn, method = method, ...)
  if (any(!is.finite(res$par)) || !is.finite(res$value))
    stop("Non-finite optimiser result")
  res
}

link_fns <- list(identity = function(z) z, softplus = softplus)

make_dist_registry <- function() {
  reg <- list(
    norm = list(
      param_names = c("mu", "sigma"),
      link_vector = c("identity", "softplus"),
      logpdf = function(x, mu, sigma) dnorm(x, mean = mu, sd = sigma, log = TRUE),
      invcdf = qnorm
    ),
    exp = list(
      param_names = "rate",
      link_vector = "softplus",
      logpdf = function(x, rate) dexp(x, rate = rate, log = TRUE),
      invcdf = qexp
    ),
    gamma = list(
      param_names = c("shape", "rate"),
      link_vector = c("softplus", "softplus"),
      logpdf = function(x, shape, rate) dgamma(x, shape = shape, rate = rate, log = TRUE),
      invcdf = qgamma
    ),
    weibull = list(
      param_names = c("shape", "scale"),
      link_vector = c("softplus", "softplus"),
      logpdf = function(x, shape, scale) dweibull(x, shape = shape, scale = scale, log = TRUE),
      invcdf = qweibull
    ),
    lnorm = list(
      param_names = c("meanlog", "sdlog"),
      link_vector = c("identity", "softplus"),
      logpdf = function(x, meanlog, sdlog) dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE),
      invcdf = qlnorm
    ),
    pois = list(
      param_names = "lambda",
      link_vector = "softplus",
      logpdf = function(x, lambda) dpois(x, lambda = lambda, log = TRUE),
      invcdf = qpois
    ),
    beta = list(
      param_names = c("shape1", "shape2"),
      link_vector = c("softplus", "softplus"),
      logpdf = function(x, shape1, shape2) dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE),
      invcdf = qbeta
    ),
    logis = list(
      param_names = c("location", "scale"),
      link_vector = c("identity", "softplus"),
      logpdf = function(x, location, scale) dlogis(x, location = location, scale = scale, log = TRUE),
      invcdf = qlogis
    )
  )
  for (n in names(reg))
    reg[[n]]$param_count <- length(reg[[n]]$param_names) * 2
  reg
}

compute_params <- function(theta, offset, family) {
  npar <- length(family$param_names)
  params <- vector("list", npar)
  for (i in seq_len(npar)) {
    b0 <- theta[2 * i - 1]
    b1 <- theta[2 * i]
    lin <- b0 + b1 * offset
    link_id <- family$link_vector[i]
    params[[i]] <- link_fns[[link_id]](lin)
  }
  params
}

make_nll <- function(family_name, x_prev, x_curr, registry = dist_registry) {
  fam <- registry[[family_name]]
  offset <- if (is.null(x_prev)) rep(0, length(x_curr)) else x_prev
  safe_support(x_curr, family_name)
  function(theta) {
    pars <- compute_params(theta, offset, fam)
    logpdf <- do.call(fam$logpdf, c(list(x_curr), pars))
    -sum(logpdf)
  }
}

nll_fun_from_cfg <- function(k, cfg, xs, Xprev, registry = dist_registry) {
  dname <- cfg[[k]]$distr
  offset <- if (k > 1) Xprev[, k - 1] else rep(0, length(xs))
  make_nll(dname, offset, xs, registry)
}

eval_ll_from_cfg <- function(k, pars, X, cfg, registry = dist_registry) {
  xs    <- X[, k]
  Xprev <- if (k > 1) X[, k - 1] else NULL
  dname <- cfg[[k]]$distr
  fam   <- registry[[dname]]
  offset <- if (k > 1) X[, k - 1] else rep(0, length(xs))
  pars_list <- compute_params(pars, offset, fam)
  do.call(fam$logpdf, c(list(xs), pars_list))
}

fit_param <- function(X_pi_train, X_pi_test, config, registry = dist_registry) {
  param_len <- sapply(config, function(cf) registry[[cf$distr]]$param_count)
  init_vals <- lapply(param_len, function(n) rep(0, n))
  param_est <- vector("list", K)
  for (k in seq_len(K)) {
    xs    <- X_pi_train[, k]
    Xprev <- if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE] else NULL
    nll   <- nll_fun_from_cfg(k, config, xs, Xprev, registry)
    fit   <- safe_optim(init_vals[[k]], nll)
    param_est[[k]] <- fit$par

    pll <- sum(eval_ll_from_cfg(k, param_est[[k]], X_pi_train, config, registry))
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
    eval_ll_from_cfg(k, param_est[[k]], X_pi_test, config, registry)
  )

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

summarise_fit <- function(param_est, X_test, ll_delta_df, cfg = config, registry = dist_registry) {
  K <- length(param_est)
  mean_param_test <- sapply(seq_len(K), function(k) {
    offset <- if (k > 1) X_test[, k - 1] else rep(0, nrow(X_test))
    pars  <- param_est[[k]]
    dname <- cfg[[k]]$distr
    fam   <- registry[[dname]]
    params <- compute_params(pars, offset, fam)
    mean(params[[1]])
  })
  mle_param <- sapply(param_est, function(p) p[1])
  ll_delta_df$mean_param_test <- round(mean_param_test, 3)
  ll_delta_df$mle_param <- round(mle_param, 3)
  ll_delta_df
}

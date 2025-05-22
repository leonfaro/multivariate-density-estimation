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
  reg
}

compute_distribution_parameters <- function(theta_vector, X_prev_matrix,
                                            family_spec, N_observations) {
  output_params_list <- list()
  current_theta_idx_start <- 1
  num_betas_per_dist_param <- if (is.null(X_prev_matrix) ||
    ncol(as.matrix(X_prev_matrix)) == 0) 1 else ncol(X_prev_matrix) + 1

  intercept_column <- rep(1.0, N_observations)
  if (is.null(X_prev_matrix) || ncol(as.matrix(X_prev_matrix)) == 0) {
    D <- matrix(intercept_column, nrow = N_observations)
  } else {
    D <- cbind(intercept_column, X_prev_matrix)
  }

  for (j in seq_along(family_spec$param_names)) {
    current_param_name <- family_spec$param_names[j]
    current_link_id <- family_spec$link_vector[j]
    selected_link_function <- link_fns[[current_link_id]]
    current_theta_idx_end <- current_theta_idx_start +
      num_betas_per_dist_param - 1
    beta_subset <- theta_vector[current_theta_idx_start:current_theta_idx_end]
    linear_predictor_values <- as.numeric(D %*% beta_subset)
    transformed_values <- selected_link_function(linear_predictor_values)
    output_params_list[[current_param_name]] <- transformed_values
    current_theta_idx_start <- current_theta_idx_end + 1
  }
  output_params_list
}

make_nll <- function(family_name, x_prev, x_curr, registry = dist_registry) {
  fam <- registry[[family_name]]
  N_obs <- length(x_curr)
  if (!is.null(x_prev)) {
    if (is.null(dim(x_prev))) {
      if (length(x_prev) != N_obs)
        stop("x_prev dimension mismatch")
      x_prev <- matrix(x_prev, ncol = 1)
    } else if (nrow(x_prev) != N_obs) {
      stop("x_prev dimension mismatch")
    }
  }
  safe_support(x_curr, family_name)
  function(theta) {
    pars <- compute_distribution_parameters(theta, x_prev, fam, N_obs)
    logpdf <- do.call(fam$logpdf, c(list(x_curr), pars))
    -sum(logpdf)
  }
}

nll_fun_from_cfg <- function(k, cfg, xs, Xprev, registry = dist_registry) {
  dname <- cfg[[k]]$distr
  X_prev <- if (k > 1) Xprev[, 1:(k - 1), drop = FALSE] else NULL
  make_nll(dname, X_prev, xs, registry)
}

eval_ll_from_cfg <- function(k, pars, X, cfg, registry = dist_registry) {
  xs    <- X[, k]
  X_prev <- if (k > 1) X[, 1:(k - 1), drop = FALSE] else NULL
  dname <- cfg[[k]]$distr
  fam   <- registry[[dname]]
  N_obs <- length(xs)
  pars_list <- compute_distribution_parameters(pars, X_prev, fam, N_obs)
  do.call(fam$logpdf, c(list(xs), pars_list))
}

fit_param <- function(X_pi_train, X_pi_test, config, registry = dist_registry) {
  K <- length(config)
  param_len <- sapply(seq_len(K), function(k) {
    fam <- registry[[config[[k]]$distr]]
    n_prev <- if (k > 1) k - 1 else 0
    (n_prev + 1) * length(fam$param_names)
  })
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
    X_prev <- if (k > 1) X_test[, 1:(k - 1), drop = FALSE] else NULL
    pars  <- param_est[[k]]
    dname <- cfg[[k]]$distr
    fam   <- registry[[dname]]
    params <- compute_distribution_parameters(pars, X_prev, fam, nrow(X_test))
    mean(params[[1]])
  })
  mle_param <- sapply(param_est, function(p) p[1])
  ll_delta_df$mean_param_test <- round(mean_param_test, 3)
  ll_delta_df$mle_param <- round(mle_param, 3)
  ll_delta_df
}

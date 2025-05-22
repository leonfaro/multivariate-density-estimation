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

make_generalized_nll <- function(family_name_str, X_prev_data_matrix,
                                 x_curr_vector, registry = dist_registry) {
  family_spec <- registry[[family_name_str]]
  safe_support(x_curr_vector, family_name_str)
  function(theta) {
    computed_params <- compute_distribution_parameters(
      theta,
      X_prev_data_matrix,
      family_spec,
      length(x_curr_vector)
    )
    logpdf_function <- family_spec$logpdf
    log_pdf_values <- do.call(logpdf_function,
                              c(list(x_curr_vector), computed_params))
    if (any(!is.finite(log_pdf_values)))
      return(Inf)
    -sum(log_pdf_values)
  }
}

fit_param <- function(X_pi_train, X_pi_test, config, registry = dist_registry) {
  K <- length(config)

  init_vals_list <- list()
  for (k_loop_var in seq_len(K)) {
    family_spec_init <- registry[[config[[k_loop_var]]$distr]]
    num_dist_params_init <- length(family_spec_init$param_names)
    num_betas_per_dist_param_init_k <-
      if (k_loop_var == 1) 1 else (k_loop_var - 1) + 1
    dynamic_theta_len_init_k <-
      num_dist_params_init * num_betas_per_dist_param_init_k
    init_vals_list[[k_loop_var]] <- rep(0, dynamic_theta_len_init_k)
  }
  init_vals <- init_vals_list

  param_est <- vector("list", K)
  for (k in seq_len(K)) {
    current_X_pi_train_k <- X_pi_train[, k]
    current_X_prev_train_matrix <-
      if (k == 1) NULL else X_pi_train[, 1:(k - 1), drop = FALSE]
    dist_name_k <- config[[k]]$distr
    nll <- make_generalized_nll(dist_name_k,
                                current_X_prev_train_matrix,
                                current_X_pi_train_k,
                                registry)
    fit <- safe_optim(init_vals[[k]], nll)
    param_est[[k]] <- fit$par

    estimated_theta_k <- param_est[[k]]
    family_spec_k <- registry[[config[[k]]$distr]]
    computed_params_train <- compute_distribution_parameters(
      estimated_theta_k,
      current_X_prev_train_matrix,
      family_spec_k,
      length(current_X_pi_train_k)
    )
    logpdf_values_train <- do.call(family_spec_k$logpdf,
                                   c(list(x = current_X_pi_train_k),
                                     computed_params_train))
    pll <- sum(logpdf_values_train)

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

  param_ll_mat_test <- sapply(seq_len(K), function(k_sapply) {
    current_X_pi_test_k <- X_pi_test[, k_sapply]
    current_X_prev_test_matrix <- if (k_sapply == 1)
      NULL else X_pi_test[, 1:(k_sapply - 1), drop = FALSE]
    family_spec_test_k <- registry[[config[[k_sapply]]$distr]]
    estimated_theta_k_from_train <- param_est[[k_sapply]]
    computed_params_test <- compute_distribution_parameters(
      estimated_theta_k_from_train,
      current_X_prev_test_matrix,
      family_spec_test_k,
      length(current_X_pi_test_k)
    )
    logpdf_values_test_k <- do.call(family_spec_test_k$logpdf,
                                    c(list(x = current_X_pi_test_k),
                                      computed_params_test))
    logpdf_values_test_k
  })

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
  mean_param_test <- sapply(seq_len(K), function(k_sapply) {
    current_X_prev_test_matrix_summary <- if (k_sapply == 1)
      NULL else X_test[, 1:(k_sapply - 1), drop = FALSE]
    num_obs_test <- nrow(X_test)
    family_spec_summary_k <- registry[[cfg[[k_sapply]]$distr]]
    estimated_theta_k_summary <- param_est[[k_sapply]]
    computed_dist_params_summary <- compute_distribution_parameters(
      estimated_theta_k_summary,
      current_X_prev_test_matrix_summary,
      family_spec_summary_k,
      num_obs_test
    )
    mean_val <- mean(computed_dist_params_summary[[1]])
    mean_val
  })
  mle_param <- sapply(param_est, function(p) p[1])
  ll_delta_df$mean_param_test <- round(mean_param_test, 3)
  ll_delta_df$mle_param <- round(mle_param, 3)
  ll_delta_df
}

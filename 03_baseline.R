# Input: X_train, X_test, cfg, dist_registry
# Output: fit_joint_param(), summary_table()
# Input: initial vector par, objective fn, method
# Output: optim result
safe_optim <- function(par, fn, method = "BFGS", ...) {
  res <- optim(par, fn, method = method, ...)
  if (any(!is.finite(res$par)) || !is.finite(res$value))
    stop("optim failed")
  res
}

# Input: positive y
# Output: inverse softplus
softplus_inv <- function(y) {
  ifelse(y > 20, y, log(exp(y) - 1))
}

# Input: cfg list, registry
# Output: list with parameter names and counts
parse_param_spec <- function(cfg, registry = dist_registry) {
  out <- list()
  param_names <- character()
  num_params <- integer(length(cfg))
  dims <- vapply(cfg, function(c) as.integer(registry[[c$distr]]$dim), integer(1))
  cum_dims <- cumsum(dims)
  for (j in seq_along(cfg)) {
    fam <- registry[[cfg[[j]]$distr]]
    p <- if (j == 1) 0 else cum_dims[j - 1]
    n_par <- length(fam$param_names)
    num_params[j] <- n_par * (p + 1)
    for (nm in fam$param_names) {
      labels <- c("intercept", if (p > 0) paste0("X", seq_len(p), "_slope"))
      dim_idx <- seq_len(dims[j]) + if (j == 1) 0 else cum_dims[j - 1]
      param_names <- c(param_names,
                       paste0("dim", rep(dim_idx, each = length(labels)), "_",
                              nm, "_", rep(labels, times = dims[j])))
    }
  }
  out$num_params <- num_params
  out$param_names <- param_names
  out$dims <- dims
  out$cum_dims <- cum_dims
  out
}

# Input: parameter vector theta, past data matrix, distribution spec, N_obs
# Output: list of transformed parameters per observation
compute_distribution_parameters <- function(theta, X_prev, family_spec, N_obs) {
  p <- if (is.null(X_prev)) 0 else ncol(X_prev)
  X_mat <- if (p == 0) matrix(0, nrow = N_obs, ncol = 0) else as.matrix(X_prev)
  out <- vector("list", length(family_spec$param_names))
  idx <- 1
  for (j in seq_along(out)) {
    betas <- theta[idx:(idx + p)]
    eta <- betas[1]
    if (p > 0)
      eta <- eta + X_mat %*% betas[-1]
    val <- link_fns[[family_spec$link_vector[j]]](eta)
    out[[j]] <- as.vector(val)
    idx <- idx + p + 1
  }
  names(out) <- family_spec$param_names
  out
}

# Input: distribution name, past data, x_vec
# Output: nll objective function
make_generalized_nll <- function(family_name, X_prev, x_vec,
                                 registry = dist_registry) {
  fam <- registry[[family_name]]
  N <- if (is.matrix(x_vec)) nrow(x_vec) else length(x_vec)
  function(theta) {
    pars <- compute_distribution_parameters(theta, X_prev, fam, N)
    ll <- do.call(fam$logpdf, c(list(x_vec), pars))
    if (any(!is.finite(ll)))
      return(Inf)
    -sum(ll)
  }
}

# Input: training matrix, test matrix, cfg list, registry
# Output: parameter estimates and log-likelihood matrices
fit_joint_param <- function(X_train, X_test, cfg, registry = dist_registry) {
  dims <- vapply(cfg, function(c) as.integer(registry[[c$distr]]$dim), integer(1))
  cum_dims <- cumsum(dims)
  K <- sum(dims)
  param_est <- vector("list", length(cfg))
  true_ll_mat <- matrix(0, nrow = nrow(X_test), ncol = K)
  joint_ll_mat <- matrix(0, nrow = nrow(X_test), ncol = K)

  for (j in seq_along(cfg)) {
    start_idx <- if (j == 1) 1 else cum_dims[j - 1] + 1
    end_idx <- cum_dims[j]
    x_tr <- X_train[, start_idx:end_idx, drop = FALSE]
    x_prev_tr <- if (start_idx > 1) X_train[, 1:(start_idx - 1), drop = FALSE] else NULL
    nll <- make_generalized_nll(cfg[[j]]$distr, x_prev_tr, x_tr, registry)
    fam <- registry[[cfg[[j]]$distr]]
    p <- if (j == 1) 0 else cum_dims[j - 1]
    init <- numeric(length(fam$param_names) * (p + 1))
    for (m in seq_along(fam$param_names)) {
      if (fam$link_vector[m] == "softplus")
        init[(m - 1) * (p + 1) + 1] <- softplus_inv(1)
    }
    param_est[[j]] <- safe_optim(init, nll)$par

    x_te <- X_test[, start_idx:end_idx, drop = FALSE]
    x_prev_te <- if (start_idx > 1) X_test[, 1:(start_idx - 1), drop = FALSE] else NULL
    pars <- compute_distribution_parameters(param_est[[j]], x_prev_te, fam,
                                            nrow(x_te))
    joint_block <- do.call(fam$logpdf, c(list(x_te), pars))
    joint_ll_mat[, start_idx:end_idx] <- joint_block
    for (kk in seq_len(dims[j])) {
      k_global <- start_idx + kk - 1
      true_ll_mat[, k_global] <- pdf_k(k_global, x_te[, kk],
                                       if (start_idx > 1) x_prev_te else numeric(0),
                                       cfg, log = TRUE)
    }
  }
  ll_df <- data.frame(
    dim = seq_len(K),
    distr = rep(sapply(cfg, `[[`, "distr"), times = dims),
    ll_true = colMeans(true_ll_mat),
    ll_joint = colMeans(joint_ll_mat)
  )
  ll_df$delta_joint <- ll_df$ll_true - ll_df$ll_joint
  ll_df[, 3:5] <- round(ll_df[, 3:5], 6)
  list(
    param_est = param_est,
    ll_delta_df_test = ll_df,
    true_ll_mat_test = true_ll_mat,
    joint_ll_mat_test = joint_ll_mat
  )
}

# Input: training data, cfg list, estimates, average log-liks, registry
# Output: data.frame summarizing parameters and likelihoods
summary_table <- function(X_train, cfg, theta_hat, LL_true_avg, LL_joint_avg,
                          registry = dist_registry) {
  thetas <- if (!is.null(theta_hat$param_est)) theta_hat$param_est else theta_hat
  dims <- vapply(cfg, function(c) as.integer(registry[[c$distr]]$dim), integer(1))
  cum_dims <- cumsum(dims)
  K <- sum(dims)
  out <- data.frame(
    dim = seq_len(K),
    distr = rep(sapply(cfg, `[[`, "distr"), times = dims),
    ll_true_avg = LL_true_avg,
    ll_joint_avg = LL_joint_avg,
    delta_joint = LL_true_avg - LL_joint_avg,
    stringsAsFactors = FALSE
  )
  p1 <- numeric(K)
  p2 <- character(K)
  m1 <- numeric(K)
  m2 <- character(K)
  for (j in seq_along(cfg)) {
    start_idx <- if (j == 1) 1 else cum_dims[j - 1] + 1
    fam <- registry[[cfg[[j]]$distr]]
    X_prev <- if (start_idx > 1) X_train[, 1:(start_idx - 1), drop = FALSE] else NULL
    pars <- compute_distribution_parameters(thetas[[j]], X_prev, fam, nrow(X_train))
    X0 <- if (start_idx > 1) matrix(0, nrow = 1, ncol = start_idx - 1) else NULL
    pars0 <- compute_distribution_parameters(thetas[[j]], X0, fam, 1)
    for (kk in seq_len(dims[j])) {
      k_gl <- start_idx + kk - 1
      p1[k_gl] <- mean(pars[[1]])
      p2[k_gl] <- if (length(pars) >= 2)
        sprintf("%.6f", mean(pars[[2]])) else "none"
      m1[k_gl] <- pars0[[1]][1]
      m2[k_gl] <- if (length(pars0) >= 2)
        sprintf("%.6f", pars0[[2]][1]) else "none"
    }
  }
  out$true_param1 <- round(p1, 6)
  out$mean_param2 <- p2
  out$mle_base1 <- round(m1, 6)
  out$mle_base2 <- m2
  out
}



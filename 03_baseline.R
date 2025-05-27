safe_optim <- function(par, fn, method = "BFGS", ...) {
  res <- optim(par, fn, method = method, ...)
  if (any(!is.finite(res$par)) || !is.finite(res$value))
    stop("optim failed")
  res
}

parse_param_spec <- function(cfg, registry = dist_registry) {
  out <- list()
  param_names <- character()
  num_params <- integer(length(cfg))
  for (k in seq_along(cfg)) {
    fam <- registry[[cfg[[k]]$distr]]
    p <- if (k == 1) 0 else k - 1
    n_par <- length(fam$param_names)
    num_params[k] <- n_par * (p + 1)
    for (nm in fam$param_names) {
      labels <- c("intercept", if (p > 0) paste0("X", seq_len(p), "_slope"))
      param_names <- c(param_names,
                       paste0("dim", k, "_", nm, "_", labels))
    }
  }
  param_names_global <<- param_names
  out$num_params <- num_params
  out
}

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

fit_joint_param <- function(X_train, X_test, cfg, registry = dist_registry) {
  K <- ncol(X_train)
  param_est <- vector("list", K)
  true_ll_mat <- matrix(0, nrow = nrow(X_test), ncol = K)
  joint_ll_mat <- true_ll_mat
  for (k in seq_len(K)) {
    x_tr <- X_train[, k]
    x_prev_tr <- if (k > 1) X_train[, 1:(k - 1), drop = FALSE] else NULL
    nll <- make_generalized_nll(cfg[[k]]$distr, x_prev_tr, x_tr, registry)
    fam <- registry[[cfg[[k]]$distr]]
    p <- if (k == 1) 0 else k - 1
    init <- rep(0, length(fam$param_names) * (p + 1))
    param_est[[k]] <- safe_optim(init, nll)$par

    x_te <- X_test[, k]
    x_prev_te <- if (k > 1) X_test[, 1:(k - 1), drop = FALSE] else NULL
    fam <- registry[[cfg[[k]]$distr]]
    pars <- compute_distribution_parameters(param_est[[k]], x_prev_te, fam,
                                            length(x_te))
    joint_ll_mat[, k] <- do.call(fam$logpdf, c(list(x_te), pars))
    true_ll_mat[, k] <- pdf_k(k, x_te, if (k > 1) x_prev_te else numeric(0),
                              cfg, log = TRUE)
  }
  ll_df <- data.frame(
    dim = seq_len(K),
    distr = sapply(cfg, `[[`, "distr"),
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

summary_table <- function(X_train, cfg, theta_hat, LL_true_avg, LL_joint_avg,
                          registry = dist_registry) {
  thetas <- if (!is.null(theta_hat$param_est)) theta_hat$param_est else theta_hat
  K <- length(thetas)
  out <- data.frame(
    dim = seq_len(K),
    distr = sapply(cfg, `[[`, "distr"),
    ll_true_avg = LL_true_avg,
    ll_joint_avg = LL_joint_avg,
    delta_joint = LL_true_avg - LL_joint_avg,
    stringsAsFactors = FALSE
  )
  p1 <- numeric(K)
  p2 <- character(K)
  m1 <- numeric(K)
  m2 <- character(K)
  for (k in seq_len(K)) {
    X_prev <- if (k > 1) X_train[, 1:(k - 1), drop = FALSE] else NULL
    fam <- registry[[out$distr[k]]]
    pars <- compute_distribution_parameters(thetas[[k]], X_prev, fam,
                                            nrow(X_train))
    p1[k] <- mean(pars[[1]])
    p2[k] <- if (length(pars) >= 2)
      sprintf("%.6f", mean(pars[[2]])) else "none"

    X0 <- if (k > 1) matrix(0, nrow = 1, ncol = k - 1) else NULL
    pars0 <- compute_distribution_parameters(thetas[[k]], X0, fam, 1)
    m1[k] <- pars0[[1]][1]
    m2[k] <- if (length(pars0) >= 2)
      sprintf("%.6f", pars0[[2]][1]) else "none"
  }
  out$true_param1 <- round(p1, 6)
  out$mean_param2 <- p2
  out$mle_base1 <- round(m1, 6)
  out$mle_base2 <- m2
  out
}



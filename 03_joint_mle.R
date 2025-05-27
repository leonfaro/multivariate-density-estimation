#' Joint MLE for triangular models
#'
#' Works for arbitrary K. To add a new dimension:
#'   1. Append a list to `config` with `$distr` and (optional) `$parm`.
#'   2. Ensure `slice_pars()` recognises any new global parameters.
#' No other code changes are required.

param_names_global <- character()

parse_param_spec <- function(config, registry) {
  fam_seen <- character()
  names_out <- character()
  for (k in seq_along(config)) {
    dname <- config[[k]]$distr
    if (dname %in% fam_seen) next
    fam_seen <- c(fam_seen, dname)
    fam_names <- switch(dname,
      norm = c("mu", "log_sigma"),
      exp = c("alpha", "beta"),
      beta = c("alpha1", "beta1", "alpha2", "beta2"),
      gamma = c("shape_alpha", "shape_beta", "rate_alpha", "rate_beta"),
      weibull = c("wshape_alpha", "wshape_beta", "wscale_alpha", "wscale_beta"),
      stop(sprintf("Unsupported distribution: %s", dname))
    )
    names_out <- c(names_out, fam_names)
  }
  assign("param_names_global", names_out, envir = .GlobalEnv)
  stopifnot(all(unique(param_names_global) == param_names_global))
  invisible(param_names_global)
}

theta2list <- function(theta) {
  stopifnot(length(theta) == length(param_names_global))
  as.list(stats::setNames(theta, param_names_global))
}

slice_pars <- function(k, theta, X_prev) {
  dname <- config[[k]]$distr
  x_last <- if (length(X_prev)) X_prev[, ncol(X_prev)] else 0
  switch(dname,
    norm = list(mu = theta$mu, sigma = exp(theta$log_sigma)),
    exp = list(rate = exp(theta$alpha + theta$beta * x_last)),
    beta = list(
      shape1 = exp(theta$alpha1 + theta$beta1 * x_last),
      shape2 = exp(theta$alpha2 + theta$beta2 * x_last)
    ),
    gamma = list(
      shape = exp(theta$shape_alpha + theta$shape_beta * x_last),
      rate = exp(theta$rate_alpha + theta$rate_beta * x_last)
    ),
    weibull = list(
      shape = exp(theta$wshape_alpha + theta$wshape_beta * x_last),
      scale = exp(theta$wscale_alpha + theta$wscale_beta * x_last)
    ),
    stop(sprintf("slice_pars not implemented for %s", dname))
  )
}

safe_optim <- function(par, fn, method = "BFGS", ...) {
  res <- optim(par, fn, method = method, ...)
  if (any(!is.finite(res$par)) || !is.finite(res$value))
    stop("Non-finite optimiser result")
  res
}

joint_loglik <- function(theta, X, config, registry) {
  theta_named <- theta2list(theta)
  n <- nrow(X)
  K <- ncol(X)
  ll <- 0
  for (i in seq_len(n)) {
    for (k in seq_len(K)) {
      X_prev <- if (k > 1) X[i, 1:(k - 1), drop = FALSE] else matrix(0, nrow = 1, ncol = 0)
      pars <- slice_pars(k, theta_named, X_prev)
      lp <- do.call(registry[[config[[k]]$distr]]$logpdf,
                    c(list(x = X[i, k]), pars))
      if (!is.finite(lp))
        return(Inf)
      ll <- ll + lp
    }
  }
  -ll
}

fit_joint_param <- function(X_train, X_test, config, registry = dist_registry) {
  parse_param_spec(config, registry)
  init <- rep(0, length(param_names_global))
  nll <- function(th) joint_loglik(th, X_train, config, registry)
  opt <- safe_optim(init, nll)
  theta_hat <- opt$par
  theta_list <- theta2list(theta_hat)

  true_ll_mat_test <- sapply(seq_len(ncol(X_test)), function(k)
    pdf_k(k, X_test[, k], if (k > 1) X_test[, 1:(k - 1), drop = FALSE] else numeric(0),
          config, log = TRUE))
  joint_ll_mat_test <- matrix(NA_real_, nrow = nrow(X_test), ncol = ncol(X_test))
  for (i in seq_len(nrow(X_test))) {
    for (k in seq_len(ncol(X_test))) {
      X_prev <- if (k > 1) X_test[i, 1:(k - 1), drop = FALSE] else matrix(0, nrow = 1, ncol = 0)
      pars <- slice_pars(k, theta_list, X_prev)
      joint_ll_mat_test[i, k] <- do.call(registry[[config[[k]]$distr]]$logpdf,
                                         c(list(x = X_test[i, k]), pars))
    }
  }
  ll_delta_df_test <- data.frame(
    dim = seq_len(ncol(X_test)),
    distr = sapply(config, `[[`, "distr"),
    ll_true = colMeans(true_ll_mat_test),
    ll_joint = colMeans(joint_ll_mat_test)
  )
  ll_delta_df_test$delta_joint <- ll_delta_df_test$ll_true - ll_delta_df_test$ll_joint
  ll_delta_df_test[, 3:5] <- round(ll_delta_df_test[, 3:5], 6)
  list(theta_hat = theta_hat, ll_delta_df_test = ll_delta_df_test)
}


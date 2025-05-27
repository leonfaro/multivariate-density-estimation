safe_optim <- function(par, fn, method = "BFGS", ...) {
  res <- optim(par, fn, method = method, ...)
  if (any(!is.finite(res$par)) || !is.finite(res$value))
    stop("Non-finite optimiser result")
  res
}

# [#1] Global theta specification
parse_param_spec <- function(config, registry) {
  fam_seen <- character()
  param_names <- character()
  for (k in seq_len(length(config))) {
    dname <- config[[k]]$distr
    if (dname %in% fam_seen) next
    fam_seen <- c(fam_seen, dname)
    fam_names <- switch(dname,
      norm   = c("mu", "log_sigma"),
      exp    = c("alpha", "beta"),
      beta   = c("alpha1", "beta1", "alpha2", "beta2"),
      gamma  = c("shape_alpha", "shape_beta", "rate_alpha", "rate_beta"),
      weibull = c("wshape_alpha", "wshape_beta", "wscale_alpha", "wscale_beta"),
      logis  = c("loc_alpha", "loc_beta", "scale_alpha", "scale_beta"),
      stop("parse_param_spec: unknown distr"))
    param_names <- c(param_names, fam_names)
  }
  unique(param_names)
}

param_names_global <- parse_param_spec(config, dist_registry)

# [#2] Helpers
theta2list <- function(theta) setNames(as.list(theta), param_names_global)

slice_pars <- function(k, theta_named, X_prev_row, cfg) {
  dname <- cfg[[k]]$distr
  x_last <- if (ncol(X_prev_row)) X_prev_row[, ncol(X_prev_row)] else 0
  switch(dname,
         norm   = list(mu = theta_named$mu,
                       sigma = exp(theta_named$log_sigma)),
         exp    = list(rate = exp(theta_named$alpha +
                                  theta_named$beta * x_last)),
         beta   = list(shape1 = exp(theta_named$alpha1 +
                                    theta_named$beta1 * x_last),
                       shape2 = exp(theta_named$alpha2 +
                                    theta_named$beta2 * x_last)),
         stop("slice_pars: unknown distr"))
}

# [#3] joint negative log-likelihood
joint_nll <- function(theta_vec, X, cfg, registry) {
  th <- theta2list(theta_vec)
  n <- nrow(X); K <- ncol(X)
  acc <- 0
  for (i in seq_len(n)) {
    for (k in seq_len(K)) {
      pars <- slice_pars(k, th, X[i, 1:(k-1), drop = FALSE], cfg)
      lp <- do.call(registry[[cfg[[k]]$distr]]$logpdf,
                   c(list(x = X[i, k]), pars))
      acc <- acc - lp
    }
  }
  if (!is.finite(acc)) acc <- Inf
  acc
}

# [#4] fit_param
fit_param <- function(X_train, X_test, config, registry = dist_registry) {
  init <- rep(0, length(param_names_global))
  nll <- function(th) joint_nll(th, X_train, config, registry)
  opt <- safe_optim(init, nll, method = "BFGS")
  theta_hat <- theta2list(opt$par)
  ll_true  <- matrix(NA_real_, nrow = nrow(X_test), ncol = ncol(X_test))
  ll_joint <- ll_true
  for (i in seq_len(nrow(X_test))) {
    for (k in seq_len(ncol(X_test))) {
      pars <- slice_pars(k, theta_hat,
                         X_test[i, 1:(k-1), drop = FALSE], config)
      ll_true[i, k]  <- pdf_k(k, X_test[i, k],
                              if (k > 1) X_test[i, 1:(k-1), drop = FALSE]
                              else numeric(0),
                              config, log = TRUE)
      ll_joint[i, k] <- do.call(registry[[config[[k]]$distr]]$logpdf,
                               c(list(x = X_test[i, k]), pars))
    }
  }
  out <- data.frame(
    dim      = seq_len(ncol(X_test)),
    distr    = sapply(config, `[[`, "distr"),
    ll_true  = colMeans(ll_true),
    ll_joint = colMeans(ll_joint)
  )
  out$delta_joint <- out$ll_true - out$ll_joint
  list(theta_hat = theta_hat,
       ll_delta_df_test = out)
}

# [#5] simple self-test
if (Sys.getenv("RUN_MLE_SELFTEST") == "1") {
  stopifnot(max(abs(fit_param(X_sim, X_sim, cfg)$ll_delta_df_test$delta_joint)) < 1e-2)
}


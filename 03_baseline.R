
# **Parametric Optimization – Detail-Workflow**
# *(R-Codex-ready: minimal prose, maximal formal spec)*
#
# | **S**  | **Function Name ** | **Input (type / shape)**                  | **Output**                                                                       | **Core Logic / Formula (precise, no prose)**                                                                                                                                         |
# | ------ | ---------------------------------------------- | ----------------------------------------- | -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
# | **-1** | **`cfg_init()`**                               | —                                         | `cfg` : list [K]                                                                | for each *k*: list(distr, n_par, parm=NULL)                                                                                                                                         |
# | **-0** | **`sample_pi_df()`**                           | `N, cfg`                                  | `X ∈ ℝ^{N×K}`                                                                    | draw `Z ~ N(0,I_K)` → `X ← S_inv(Z,cfg)`                                                                                                                                             |
# | **0**  | **`split_X()`**                                | `X`                                       | list `X_k, X_prev`                                                               | `X_prev ← X[,1:(k-1)]`, `X_k ← X[,k]`                                                                                                                                                |
# | **1**  | **`design_k()`**                               | `X_prev ∈ ℝ^{N×(k-1)}`                    | `D_k ∈ ℝ^{N×P_k}`                                                                | `D_k ← cbind(1, X_prev)` (*all* predecessors)                                                                                                                                        |
# | **2**  | **`link_k()`**                                 | `theta_k ∈ ℝ^{P_k·M_k}, D_k, distr_k`     | `pars_i`                                                                         | *Split θ_k into M_k blocks (per density):*<br>• normal: `(β_μ, β_σ)`<br>• gamma: `(β_shape, β_rate)`<br>`η_j ← D_k β_j`;<br>`param_j ← softplus(η_j)` if strictly >0 else identity |
# | **3**  | **`logpdf_k_vec()`**                           | `x_k, pars_i, distr_k`                    | `ℓ_{k,·} ∈ ℝ^N`                                                                  | element-wise `log f_{distr_k}(x_k; pars_i)`                                                                                                                                          |
# | **4**  | **`logLik_k()`**                               | `theta_k, X, k, cfg`                      | `ℓ_k ∈ ℝ`                                                                        | `D_k ← design_k(X_prev)`;<br>`pars ← link_k(theta_k,D_k,distr_k)`;<br>`ℓ_k ← sum( logpdf_k_vec(x_k,pars,distr_k) )`                                                                  |
# | **5**  | **`optim_k()`**                                | `k, start_k, X, cfg`                      | `thetâ_k`                                                                       | `optim(par=start_k, fn=function(θ) -logLik_k(θ,X,k,cfg), method="Nelder-Mead")`                                                                                                      |
# | **6**  | **`update_cfg()`**                             | `k, thetâ_k, cfg`                        | `cfg*`                                                                           | `cfg[[k]]$parm ← function(X_prev) { D←design_k(X_prev); link_k(thetâ_k,D,cfg[[k]]$distr) }`                                                                                         |
# | **7**  | **`fit_seq()`**                                | `X_train, cfg, starts`                    | list `(thetâ₁:K, cfg*)`                                                         | `for k in 1:K { thetâ_k←optim_k(...); cfg←update_cfg(...) }`                                                                                                                        |
# | **8**  | **`logLik_matrix()`**                          | `X, cfg`                                  | `L ∈ ℝ^{N×K}`                                                                    | loop k: row ℓ_{k,·} via steps 1–3                                                                                                                                                   |
# | **9a** | **`LL_train()`**                               | `X_train, cfg*`                           | `LL_train_avg`                                                                   | `mean( rowSums( logLik_matrix(X_train,cfg*) ) )`                                                                                                                                     |
# | **9b** | **`LL_test_param()`**                          | `X_test, cfg*`                            | `LL_param_avg`                                                                   | same formula on test batch                                                                                                                                                           |
# | **9c** | **`LL_test_true()`**                           | `X_test, cfg_true`                        | `LL_true_avg`                                                                    | use ground-truth cfg_true                                                                                                                                                           |
# | **10** | **`metrics()`**                                | `LL_true_avg, LL_param_avg`               | `Δ = (LL_true_avg − LL_param_avg)`                                               | expect Δ≈0 under perfect fit                                                                                                                                                         |
# | **11** | **`summary_table()`**                          | `(X_train ∈ ℝ^{N×K}, cfg* list[K], θ̂ list[K], LL_true_avg, LL_param_avg)` | `data.frame` columns →<br>• `dim` (k)<br>• `distr` (cfg*[[k]])<br>• `ll_true_avg` (= LL_true_avg)<br>• `ll_param_avg` (= LL_param_avg)<br>• `delta` (= ll_true_avg − ll_param_avg)<br>• `mean_param1`, `mean_param2`<br>• `mle_param1`, `mle_param2` | **For each k = 1…K**<br>1. `X_prev ← X_train[,1:(k-1)]` ; `D_k ← design_k(X_prev)`<br>2. `Pars ← link_k(θ̂_k, D_k, distr_k)` ⇒ matrix N×M_k (M_k = #param-groups)<br>3. `mean_param1 ← mean(Pars[,1])`<br>  `mean_param2 ← if (M_k≥2) mean(Pars[,2]) else "none"`<br>4. `D_ref ← c(1, rep(0,k-1))` ; `Pars_ref ← link_k(θ̂_k, matrix(D_ref, nrow=1), distr_k)`<br>5. `mle_param1 ← Pars_ref[1]`<br>  `mle_param2 ← if (M_k≥2) Pars_ref[2] else "none"`<br>6. append row with above values |

#
# ### Data-Shapes & Constants
#
# ```text
# N           : sample size per batch (train = test)
# K           : dimension of X
# P_k         : k (incl. intercept)  ⇒ P_k = k
# M_k         : #parameter-groups for density k
#               normal→2, gamma→2, beta→2, …
# theta_k     : length = P_k · M_k
# design_k    : N × P_k
# ```
#
# ### Density Blocks (built-in)
#
# | `distr`  | param-groups *(η_j)* | pdf formula                                                       |
# | -------- | --------------------- | ----------------------------------------------------------------- |
# | `normal` | μ, σ > 0              | `log ϕ((x-μ)/σ) - log σ`                                          |
# | `gamma`  | shape > 0, rate > 0   | `shape·log(rate) + (shape-1)·log x - rate·x - log Γ(shape)`       |
# | `beta`   | α > 0, β > 0          | `log Γ(α+β) - log Γ(α) - log Γ(β) + (α-1)·log x + (β-1)·log(1-x)` |
#
# ### Overall Workflow (single line)
#
# ```
# X_train,X_test ← sample_pi_df(N,cfg_true)
# thetâ, cfg*   ← fit_seq(X_train,cfg_init(),starts)
# LL_train_avg   ← LL_train(X_train,cfg*)
# LL_param_avg   ← LL_test_param(X_test,cfg*)
# LL_true_avg    ← LL_test_true(X_test,cfg_true)
# Δ              ← metrics(LL_true_avg,LL_param_avg)
# tbl            ← summary_table(cfg*,thetâ,LL_true_avg,LL_param_avg)
# ```
#
# All steps are strictly sequential in *k*, use full predecessor set $X_{1:k-1}$, require no gradients, and deliver per-dimension MLEs plus out-of-sample log-likelihood diagnostics normalized by *N*.

# ----- Implementation -----

safe_optim <- function(par, fn, method = "BFGS", ...) {
  res <- optim(par, fn, method = method, ...)
  if (any(!is.finite(res$par)) || !is.finite(res$value))
    stop("Non-finite optimiser result")
  res
}

compute_distribution_parameters <- function(theta_vector, X_prev_matrix,
                                            family_spec, N_observations) {
  if (is.null(X_prev_matrix) || ncol(as.matrix(X_prev_matrix)) == 0) {
    D <- matrix(1, nrow = N_observations)
  } else {
    D <- cbind(1, X_prev_matrix)
  }
  num_betas <- ncol(D)
  res <- list()
  idx <- 1
  for (j in seq_along(family_spec$param_names)) {
    idx_end <- idx + num_betas - 1
    beta_subset <- theta_vector[idx:idx_end]
    eta <- as.numeric(D %*% beta_subset)
    link_name <- family_spec$link_vector[j]
    res[[family_spec$param_names[j]]] <- link_fns[[link_name]](eta)
    idx <- idx_end + 1
  }
  res
}

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
      logis = c("loc_alpha", "loc_beta", "scale_alpha", "scale_beta"),
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

slice_pars <- function(k, theta, X_prev, cfg) {
  dname <- cfg[[k]]$distr
  x_last <- if (length(X_prev)) X_prev[, ncol(as.matrix(X_prev))] else 0
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
    logis = list(
      location = theta$loc_alpha + theta$loc_beta * x_last,
      scale = exp(theta$scale_alpha + theta$scale_beta * x_last)
    ),
    stop(sprintf("slice_pars not implemented for %s", dname))
  )
}

transform_theta <- function(theta, spec) {
  pars <- vector("list", length(spec$param_names))
  for (j in seq_along(spec$param_names)) {
    pars[[j]] <- link_fns[[spec$link_vector[j]]](theta[j])
  }
  names(pars) <- spec$param_names
  pars
}

nll_distribution <- function(theta, x, spec) {
  pars <- transform_theta(theta, spec)
  -sum(do.call(spec$logpdf, c(list(x = x), pars)))
}

mle_distribution <- function(x, dname, registry) {
  spec <- registry[[dname]]
  start <- rep(0, length(spec$param_names))
  res <- safe_optim(start, nll_distribution, spec = spec, x = x)
  list(theta_hat = res$par, pars_hat = transform_theta(res$par, spec))
}

fit_joint_param <- function(X_train, X_test, config, registry = dist_registry) {
  K <- ncol(X_train)
  est_list <- vector("list", K)
  true_ll_mat_test <- matrix(NA_real_, nrow = nrow(X_test), ncol = K)
  joint_ll_mat_test <- matrix(NA_real_, nrow = nrow(X_test), ncol = K)
  for (k in seq_len(K)) {
    dname <- config[[k]]$distr
    est_list[[k]] <- mle_distribution(X_train[, k], dname, registry)
    true_ll_mat_test[, k] <- pdf_k(k, X_test[, k],
      if (k > 1) X_test[, 1:(k - 1), drop = FALSE] else numeric(0),
      config, log = TRUE)
    joint_ll_mat_test[, k] <- do.call(registry[[dname]]$logpdf,
      c(list(x = X_test[, k]), est_list[[k]]$pars_hat))
  }
  ll_delta_df_test <- data.frame(
    dim = seq_len(K),
    distr = sapply(config, `[[`, "distr"),
    ll_true = colMeans(true_ll_mat_test),
    ll_joint = colMeans(joint_ll_mat_test)
  )
  ll_delta_df_test$delta_joint <- ll_delta_df_test$ll_true - ll_delta_df_test$ll_joint
  ll_delta_df_test[, 3:5] <- round(ll_delta_df_test[, 3:5], 6)
  list(theta_hat = lapply(est_list, `[[`, "theta_hat"),
       pars_hat = lapply(est_list, `[[`, "pars_hat"),
       ll_delta_df_test = ll_delta_df_test,
       true_ll_mat_test = true_ll_mat_test,
       joint_ll_mat_test = joint_ll_mat_test)
}

summary_table <- function(X_train, cfg, fit_res,
                          LL_true_avg, LL_joint_avg,
                          registry = dist_registry) {
  K <- length(cfg)
  out <- data.frame(
    dim = seq_len(K),
    distr = sapply(cfg, `[[`, "distr"),
    ll_true_avg = LL_true_avg,
    ll_joint_avg = LL_joint_avg,
    delta_joint = LL_true_avg - LL_joint_avg,
    stringsAsFactors = FALSE
  )
  mean_p1 <- numeric(K)
  mean_p2 <- character(K)
  mle_p1 <- numeric(K)
  mle_p2 <- character(K)
  for (k in seq_len(K)) {
    X_prev <- if (k > 1) X_train[, 1:(k - 1), drop = FALSE] else matrix(0, nrow = nrow(X_train), ncol = 0)
    true_pars <- get_pars(k, X_prev, cfg)
    if (length(true_pars) == 0) {
      mean_p1[k] <- NA_real_
      mean_p2[k] <- "none"
    } else {
      mean_p1[k] <- mean(true_pars[[1]])
      if (length(true_pars) >= 2) {
        mean_p2[k] <- sprintf("%.6f", mean(true_pars[[2]]))
      } else {
        mean_p2[k] <- "none"
      }
    }
    est_pars <- fit_res$pars_hat[[k]]
    mle_p1[k] <- est_pars[[1]]
    if (length(est_pars) >= 2) {
      mle_p2[k] <- sprintf("%.6f", est_pars[[2]])
    } else {
      mle_p2[k] <- "none"
    }
  }
  out$true_param1 <- round(mean_p1, 6)
  out$mean_param2 <- mean_p2
  out$mle_base1 <- round(mle_p1, 6)
  out$mle_base2 <- mle_p2
  out
}



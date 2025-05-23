
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

make_generalized_nll <- function(family_name_str, X_prev_data_matrix,
                                 x_curr_vector, registry = dist_registry) {
  family_spec <- registry[[family_name_str]]
  safe_support(x_curr_vector, family_name_str)
  function(theta) {
    pars <- compute_distribution_parameters(theta, X_prev_data_matrix,
                                            family_spec, length(x_curr_vector))
    logpdf_vals <- do.call(family_spec$logpdf, c(list(x_curr_vector), pars))
    if (any(!is.finite(logpdf_vals)))
      return(Inf)
    -sum(logpdf_vals)
  }
}

fit_param <- function(X_pi_train, X_pi_test, config, registry = dist_registry) {
  K <- length(config)
  param_len <- sapply(seq_len(K), function(k) {
    fam <- registry[[config[[k]]$distr]]
    (if (k > 1) k else 1) * length(fam$param_names)
  })
  init_vals <- lapply(param_len, function(n) rep(0, n))
  param_est <- vector("list", K)
  for (k in seq_len(K)) {
    xs <- X_pi_train[, k]
    X_prev <- if (k > 1) X_pi_train[, 1:(k - 1), drop = FALSE] else NULL
    dname <- config[[k]]$distr
    nll <- make_generalized_nll(dname, X_prev, xs, registry)
    fit <- safe_optim(init_vals[[k]], nll)
    param_est[[k]] <- fit$par
  }
  true_ll_mat_test <- sapply(seq_len(K), function(k)
    pdf_k(k, X_pi_test[, k],
          if (k > 1) X_pi_test[, 1:(k - 1), drop = FALSE] else numeric(0),
          config, log = TRUE))
  param_ll_mat_test <- sapply(seq_len(K), function(k) {
    xs <- X_pi_test[, k]
    X_prev <- if (k > 1) X_pi_test[, 1:(k - 1), drop = FALSE] else NULL
    dname <- config[[k]]$distr
    fam <- registry[[dname]]
    pars <- compute_distribution_parameters(param_est[[k]], X_prev, fam, length(xs))
    do.call(fam$logpdf, c(list(xs), pars))
  })
  ll_delta_df_test <- data.frame(
    dim = seq_len(K),
    distribution = sapply(config, `[[`, "distr"),
    ll_true_avg = apply(true_ll_mat_test, 2, mean),
    ll_param_avg = colMeans(param_ll_mat_test)
  )
  ll_delta_df_test$delta_ll_param_avg <-
    ll_delta_df_test$ll_true_avg - ll_delta_df_test$ll_param_avg
  ll_delta_df_test[, 3:5] <- round(ll_delta_df_test[, 3:5], 6)
  list(param_est = param_est,
       ll_delta_df_test = ll_delta_df_test,
       true_ll_mat_test = true_ll_mat_test,
       param_ll_mat_test = param_ll_mat_test)
}

summary_table <- function(X_train, cfg, theta_hat,
                          LL_true_avg, LL_param_avg,
                          registry = dist_registry) {
  K <- length(theta_hat)
  out <- data.frame(
    dim = seq_len(K),
    distr = sapply(cfg, `[[`, "distr"),
    ll_true_avg = LL_true_avg,
    ll_param_avg = LL_param_avg,
    delta = LL_true_avg - LL_param_avg,
    stringsAsFactors = FALSE
  )
  mean_p1 <- numeric(K)
  mean_p2 <- character(K)
  mle_p1 <- numeric(K)
  mle_p2 <- character(K)
  for (k in seq_len(K)) {
    X_prev <- if (k > 1) X_train[, 1:(k - 1), drop = FALSE] else NULL
    fam <- registry[[out$distr[k]]]
    pars <- compute_distribution_parameters(theta_hat[[k]], X_prev,
                                            fam, nrow(X_train))
    mean_p1[k] <- mean(pars[[1]])
    if (length(pars) >= 2) {
      mean_p2[k] <- sprintf("%.6f", mean(pars[[2]]))
    } else {
      mean_p2[k] <- "none"
    }
    X0 <- if (k > 1) matrix(0, nrow = 1, ncol = k - 1) else NULL
    pars_ref <- compute_distribution_parameters(theta_hat[[k]], X0,
                                               fam, 1)
    mle_p1[k] <- pars_ref[[1]][1]
    if (length(pars_ref) >= 2) {
      mle_p2[k] <- sprintf("%.6f", pars_ref[[2]][1])
    } else {
      mle_p2[k] <- "none"
    }
  }
  out$mean_param1 <- round(mean_p1, 6)
  out$mean_param2 <- mean_p2
  out$mle_param1 <- round(mle_p1, 6)
  out$mle_param2 <- mle_p2
  out
}

save_estimated_betas <- function(param_est_list, config_list,
                                dist_registry_obj, output_dir_path) {
  if (!dir.exists(output_dir_path))
    dir.create(output_dir_path, recursive = TRUE)
  for (k in seq_along(config_list)) {
    dist_name_k <- config_list[[k]]$distr
    est_theta_k <- param_est_list[[k]]
    family_spec_k <- dist_registry_obj[[dist_name_k]]
    npar <- length(family_spec_k$param_names)
    nbeta <- if (k == 1) 1 else (k - 1) + 1
    beta_names_k <- character(length(est_theta_k))
    idx <- 1
    for (p_idx in seq_len(npar)) {
      beta_names_k[idx] <- paste0(family_spec_k$param_names[p_idx], "_intercept")
      idx <- idx + 1
      if (nbeta > 1) {
        for (j_prev in seq_len(k - 1)) {
          beta_names_k[idx] <- paste0(family_spec_k$param_names[p_idx],
                                     "_X", j_prev, "_slope")
          idx <- idx + 1
        }
      }
    }
    df <- data.frame(beta_value = est_theta_k,
                     beta_name = beta_names_k)
    csv_file <- file.path(output_dir_path,
                          paste0("dim", k, "_", dist_name_k,
                                 "_estimated_betas.csv"))
    write.csv(df, csv_file, row.names = FALSE)
  }
}

save_detailed_comparison_data <- function(data_matrix, param_ests_list,
                                          config_list, dist_registry_obj,
                                          output_dir_path,
                                          data_label_str = c("train", "test")) {
  data_label_str <- match.arg(data_label_str, c("train", "test"))
  if (!dir.exists(output_dir_path))
    dir.create(output_dir_path, recursive = TRUE)
  N_obs <- nrow(data_matrix)
  K_dims <- ncol(data_matrix)
  for (k in seq_len(K_dims)) {
    dist_name_k <- config_list[[k]]$distr
    family_spec_k <- dist_registry_obj[[dist_name_k]]
    est_theta_k <- param_ests_list[[k]]
    X_prev <- if (k == 1) NULL else data_matrix[, 1:(k - 1), drop = FALSE]
    X_k_vec <- data_matrix[, k]
    true_pars_list <- list()
    if (is.null(config_list[[k]]$parm)) {
      if (dist_name_k == "norm") {
        if ("mu" %in% family_spec_k$param_names)
          true_pars_list[["mu"]] <- rep(0, N_obs)
        if ("sigma" %in% family_spec_k$param_names)
          true_pars_list[["sigma"]] <- rep(1, N_obs)
      } else {
        for (p_name in family_spec_k$param_names)
          true_pars_list[[p_name]] <- rep(NA_real_, N_obs)
      }
      true_params_table <- as.data.frame(true_pars_list)
      if (ncol(true_params_table) > 0)
        true_params_table <- true_params_table[, family_spec_k$param_names, drop = FALSE]
    } else {
      true_map <- get_pars(k, if (k == 1) NULL else X_prev, config_list)
      true_params_table <- as.data.frame(true_map)
      true_params_table <- true_params_table[, family_spec_k$param_names, drop = FALSE]
    }
    if (ncol(true_params_table) > 0)
      colnames(true_params_table) <- paste0("true_", family_spec_k$param_names)
    fitted_params <- compute_distribution_parameters(est_theta_k, X_prev,
                                                     family_spec_k, N_obs)
    fitted_params_table <- as.data.frame(fitted_params)
    fitted_params_table <- fitted_params_table[, family_spec_k$param_names, drop = FALSE]
    colnames(fitted_params_table) <- paste0("fitted_", family_spec_k$param_names)
    true_ll <- pdf_k(k, X_k_vec, if (k == 1) numeric(0) else X_prev,
                     config_list, log = TRUE)
    fitted_ll <- do.call(family_spec_k$logpdf,
                         c(list(x = X_k_vec), fitted_params))
    output_df <- data.frame(
      X_k_value = X_k_vec,
      true_params_table,
      fitted_params_table,
      true_log_pdf = true_ll,
      fitted_log_pdf = fitted_ll,
      check.names = FALSE
    )
    csv_name <- file.path(output_dir_path,
                          paste0("dim", k, "_", dist_name_k,
                                 "_params_logpdf_", data_label_str, ".csv"))
    write.csv(output_df, csv_name, row.names = FALSE)
  }
}

run_all_diagnostics <- function(X_train, X_test, param_ests, config_list,
                                dist_registry_obj, output_dir_base) {
  save_estimated_betas(param_ests, config_list,
                       dist_registry_obj, output_dir_base)
  save_detailed_comparison_data(X_train, param_ests, config_list,
                                dist_registry_obj, output_dir_base,
                                data_label_str = "train")
  save_detailed_comparison_data(X_test, param_ests, config_list,
                                dist_registry_obj, output_dir_base,
                                data_label_str = "test")
  message(sprintf("Diagnostic CSVs for analysis saved to: %s", output_dir_base))
}

if (sys.nframe() == 0) {
  set.seed(SEED)
  data <- generate_data()
  res <- fit_param(data$train$sample$X_pi, data$test$sample$X_pi, config)
  eval_tab <- summarise_fit(res$param_est, data$test$sample$X_pi,
                            res$ll_delta_df_test, config)
  print(eval_tab)
}


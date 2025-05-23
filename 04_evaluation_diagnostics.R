# Evaluation diagnostics utilities

save_estimated_betas <- function(param_est_list, config_list,
                                dist_registry_obj, output_dir_path) {
  if (!dir.exists(output_dir_path))
    dir.create(output_dir_path, recursive = TRUE)

  for (k in seq_along(config_list)) {
    dist_name_k <- config_list[[k]]$distr
    estimated_theta_k <- param_est_list[[k]]
    family_spec_k <- dist_registry_obj[[dist_name_k]]
    num_dist_params_k <- length(family_spec_k$param_names)
    num_betas_per_dist_param_k <- if (k == 1) 1 else (k - 1) + 1
    beta_names_k <- character(length(estimated_theta_k))
    idx <- 1
    for (p_idx in seq_len(num_dist_params_k)) {
      beta_names_k[idx] <- paste0(family_spec_k$param_names[p_idx],
                                  "_intercept")
      idx <- idx + 1
      if (num_betas_per_dist_param_k > 1) {
        for (j_prev in seq_len(k - 1)) {
          beta_names_k[idx] <- paste0(family_spec_k$param_names[p_idx],
                                     "_X", j_prev, "_slope")
          idx <- idx + 1
        }
      }
    }
    df <- data.frame(beta_value = estimated_theta_k,
                     beta_name = beta_names_k)
    csv_filename <- file.path(output_dir_path,
                              paste0("dim", k, "_", dist_name_k,
                                     "_estimated_betas.csv"))
    write.csv(df, csv_filename, row.names = FALSE)
  }
}

save_detailed_comparison_data <- function(data_matrix, param_ests, config,
                                           dist_registry, output_dir,
                                           label = c("train", "test")) {
  label <- match.arg(label)
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  for (k in seq_along(config)) {
    dist_name <- config[[k]]$distr
    theta_k <- param_ests[[k]]
    family_spec <- dist_registry[[dist_name]]
    X_prev <- if (k == 1) NULL else data_matrix[, 1:(k - 1), drop = FALSE]
    X_k_obs <- data_matrix[, k]
    N_obs <- nrow(data_matrix)

    if (is.null(config[[k]]$parm)) {
      true_params_table <- data.frame(matrix(nrow = N_obs, ncol = 0))
    } else {
      true_params_map <- get_pars(k,
        if (k == 1) matrix(nrow = N_obs, ncol = 0) else X_prev,
        config)
      true_params_table <- as.data.frame(true_params_map)
      colnames(true_params_table) <-
        paste0("true_", colnames(true_params_table))
    }

    fitted_params_map <- compute_distribution_parameters(
      theta_k, X_prev, family_spec, N_obs
    )
    fitted_params_table <- as.data.frame(fitted_params_map)
    colnames(fitted_params_table) <-
      paste0("fitted_", colnames(fitted_params_table))

    true_logpdf_vec <- pdf_k(k, X_k_obs,
                             if (k == 1) numeric(0) else X_prev,
                             config, log = TRUE)
    fitted_logpdf_vec <- do.call(family_spec$logpdf,
                                 c(list(x = X_k_obs), fitted_params_map))

    output_table <- data.frame(
      X_k_obs = X_k_obs,
      true_params_table,
      fitted_params_table,
      true_logpdf = true_logpdf_vec,
      fitted_logpdf = fitted_logpdf_vec,
      check.names = FALSE
    )

    csv_name <- file.path(output_dir,
                          paste0("dim", k, "_", dist_name,
                                 "_params_logpdf_", label, ".csv"))
    write.csv(output_table, csv_name, row.names = FALSE)
  }
}


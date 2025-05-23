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
    estimated_theta_k <- param_ests_list[[k]]
    current_X_prev_matrix <- if (k == 1) NULL else data_matrix[, 1:(k - 1), drop = FALSE]
    current_X_k_vector <- data_matrix[, k]

    true_params_list_for_df <- list()
    if (is.null(config_list[[k]]$parm)) {
      if (dist_name_k == "norm") {
        if ("mu" %in% family_spec_k$param_names)
          true_params_list_for_df[["mu"]] <- rep(0.0, N_obs)
        if ("sigma" %in% family_spec_k$param_names)
          true_params_list_for_df[["sigma"]] <- rep(1.0, N_obs)
      } else {
        for (p_name in family_spec_k$param_names)
          true_params_list_for_df[[p_name]] <- rep(NA_real_, N_obs)
      }
      true_params_table <- as.data.frame(true_params_list_for_df)
      if (ncol(true_params_table) > 0)
        true_params_table <- true_params_table[, family_spec_k$param_names, drop = FALSE]
    } else {
      true_params_map_vectorized <- get_pars(k,
        if (k == 1) NULL else current_X_prev_matrix,
        config_list)
      true_params_table <- as.data.frame(true_params_map_vectorized)
      true_params_table <- true_params_table[, family_spec_k$param_names, drop = FALSE]
    }
    if (ncol(true_params_table) > 0)
      colnames(true_params_table) <- paste0("true_", family_spec_k$param_names)

    fitted_params_map_of_vectors <- compute_distribution_parameters(
      estimated_theta_k, current_X_prev_matrix, family_spec_k, N_obs
    )
    fitted_params_table <- as.data.frame(fitted_params_map_of_vectors)
    fitted_params_table <- fitted_params_table[, family_spec_k$param_names, drop = FALSE]
    colnames(fitted_params_table) <- paste0("fitted_", family_spec_k$param_names)

    true_log_pdf_vector <- pdf_k(k, current_X_k_vector,
                                 if (k == 1) numeric(0) else current_X_prev_matrix,
                                 config_list, log = TRUE)
    fitted_log_pdfs_vector <- do.call(family_spec_k$logpdf,
                                      c(list(x = current_X_k_vector),
                                        fitted_params_map_of_vectors))

    output_df <- data.frame(
      X_k_value = current_X_k_vector,
      true_params_table,
      fitted_params_table,
      true_log_pdf = true_log_pdf_vector,
      fitted_log_pdf = fitted_log_pdfs_vector,
      check.names = FALSE
    )

    csv_file_name <- file.path(output_dir_path,
                               paste0("dim", k, "_", dist_name_k,
                                      "_params_logpdf_", data_label_str, ".csv"))
    write.csv(output_df, csv_file_name, row.names = FALSE)
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

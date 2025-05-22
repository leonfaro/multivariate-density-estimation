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

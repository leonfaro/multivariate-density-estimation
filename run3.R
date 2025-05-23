

N <- 50
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] == "big") N <- 10000

config_choice <- 3
Sys.setenv(N_total = N)

source("00_setup.R")
options(width = 150)
set.seed(SEED)
data <- generate_data(N_total = N)
write_data(data)

train_df <- data$train$df
test_df  <- data$test$df
X_pi_train <- data$train$sample$X_pi
X_pi_test  <- data$test$sample$X_pi

X_train <- as.matrix(train_df[paste0("Xpi", seq_len(K))])
X_test  <- as.matrix(test_df[paste0("Xpi", seq_len(K))])

cat("Train EDA:\n")
cat("- X_pi_train Mean: ", paste(round(colMeans(X_train), 3), collapse = ", "), "\n")
cat("- X_pi_train SD:   ", paste(round(apply(X_train, 2, sd), 3), collapse = ", "), "\n")
cat("- log_det(J)_train Range: [",
    round(min(train_df$det_J), 3), ",",
    round(max(train_df$det_J), 3), "]\n\n")

cat("Test EDA:\n")
cat("- X_pi_test Mean: ", paste(round(colMeans(X_test), 3), collapse = ", "), "\n")
cat("- X_pi_test SD:   ", paste(round(apply(X_test, 2, sd), 3), collapse = ", "), "\n")
cat("- log_det(J)_test Range: [",
    round(min(test_df$det_J), 3), ",",
    round(max(test_df$det_J), 3), "]\n\n")

summary_stats <- data.frame(
  dim = seq_len(K),
  mean_train = round(colMeans(X_train), 3),
  sd_train = round(apply(X_train, 2, sd), 3),
  mean_test = round(colMeans(X_test), 3),
  sd_test = round(apply(X_test, 2, sd), 3),
  row.names = NULL
)
print(summary_stats)

param_res <- fit_param(X_pi_train, X_pi_test, config)
param_est <- param_res$param_est
tbl <- summary_table(
  X_pi_train,
  config,
  param_est,
  param_res$ll_delta_df_test$ll_true_avg,
  param_res$ll_delta_df_test$ll_param_avg
)

tbl_out <- tbl[

  , c(
    "dim", "distr", "ll_true_avg", "ll_param_avg", "delta",
    "mean_param1", "mean_param2", "mle_param1", "mle_param2"
  )
]
num_cols <- intersect(
  c("ll_true_avg", "ll_param_avg", "delta",
    "mean_param1", "mean_param2", "mle_param1", "mle_param2"),
  names(tbl_out)
)
tbl_out[num_cols] <- lapply(tbl_out[num_cols], function(x) {
  if (is.numeric(x)) sprintf("%.6f", x) else x
})
print(tbl_out, row.names = FALSE)


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

run_pipeline <- function(N_local = N) {
  Sys.setenv(N_total = N_local)
  data <- generate_data(N_total = N_local, cfg = config)
  param_res <- fit_param(data$train$sample$X_pi, data$test$sample$X_pi, config)
  tbl <- summary_table(
    data$train$sample$X_pi,
    config,
    param_res$param_est,
    param_res$ll_delta_df_test$ll_true_avg,
    param_res$ll_delta_df_test$ll_param_avg
  )
  print(tbl)
  invisible(tbl)
}

if (sys.nframe() == 0) {
  run_pipeline(N)
}

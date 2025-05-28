
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = d$X1)),
  list(distr = "beta", parm  = function(d) list(shape1 = d$X2, shape2 = 1))
)

N <- 50

source("00_setup.R")
source("01_map_definition_S.R")
source("02_sampling.R")
source("03_baseline.R")


permute_data <- function(mat, perm) mat[, perm, drop = FALSE]

run_pipeline <- function(N_local = N, cfg = config, perm = NULL) {
  Sys.setenv(N_total = N_local)
  set.seed(SEED)
  data <- generate_data(N_total = N_local, cfg = cfg)
  train_df <<- data$train$df
  test_df  <<- data$test$df
  X_pi_train <<- data$train$sample$X_pi
  X_pi_test  <<- data$test$sample$X_pi
  ll_test <<- data$test$df$ll_true
  N_test <<- nrow(X_pi_test)
  
  if (!is.null(perm)) {
    X_pi_train <<- permute_data(X_pi_train, perm)
    X_pi_test  <<- permute_data(X_pi_test, perm)
    cfg <- cfg[perm]
  }
  param_res <- fit_joint_param(X_pi_train, X_pi_test, cfg)
  ll_delta_df_test <<- param_res$ll_delta_df_test
  true_ll_mat_test <<- param_res$true_ll_mat_test
  
  source("04_forest_models.R")
  forest_res <- fit_forest(X_pi_train, X_pi_test)
  LD_hat <<- forest_res$LD_hat
  
  source("06_kernel_smoothing.R")
  ks_model <- fit_kernel(as.data.frame(X_pi_train))
  KS_hat <<- predict(ks_model, newdata = as.data.frame(X_pi_test), type = "logdensity")
  
  source("05_joint_evaluation.R", local = TRUE)
  
  target_ll <- sum(ll_test)
  forest_mismatch <<- sum(loglik_trtf) - target_ll
  kernel_mismatch <<- sum(loglik_kernel) - target_ll
  
  eval_tab <- read.csv("results/evaluation_summary.csv")
  eval_tab <- eval_tab[, !(names(eval_tab) %in%
                             c("true_param1", "mean_param2",
                               "mle_base1", "mle_base2"))]
  num_cols <- names(eval_tab)[sapply(eval_tab, is.numeric)]
  sum_row <- eval_tab[1, , drop = FALSE]
  for (col in names(sum_row)) {
    if (col %in% num_cols) {
      sum_val <- sum(abs(eval_tab[[col]]))
      if (grepl("delta", col)) {
        sum_row[[col]] <- min(sum_val, 1)
      } else {
        sum_row[[col]] <- sum_val
      }
    } else {
      sum_row[[col]] <- "sum"
    }
  }
  eval_tab <- rbind(eval_tab, sum_row)
  name <- if (is.null(perm)) "natural" else "perm"
  new_file <- file.path("results", paste0("evaluation_summary_", name, ".csv"))
  if (file.exists(new_file)) file.remove(new_file)
  file.rename("results/evaluation_summary.csv", new_file)
  print(eval_tab)
  if (is.null(perm)) {
    eval_tab_nat <<- eval_tab
  } else {
    eval_tab_perm <<- eval_tab
  }
  eval_tab <<- eval_tab
  source("dump_run5_code.R")
  invisible(eval_tab)
}


run_joint_pipeline <- function(N_local = N) {
  run_pipeline(N_local)
  joint_res <- fit_joint_param(X_pi_train, X_pi_test, config)
  eval_df$ll_joint <- joint_res$ll_delta_df_test$ll_joint
  eval_df$delta_ll_joint <- joint_res$ll_delta_df_test$delta_joint
  print(eval_df)
  invisible(eval_df)
}


eval_tab_nat  <- run_pipeline(N)
eval_tab_perm <- run_pipeline(N, perm = c(1, 3, 2))
eval_tab <- eval_tab_nat
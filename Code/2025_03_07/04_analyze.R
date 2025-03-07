# 04_analyze.R
# -------------------------------------------------------------------
# Funktionen zur Modellbewertung: Log-Likelihood, KL-Divergenz etc.
# -------------------------------------------------------------------

logLik_on_data <- function(dvals) {
  sum(log(pmax(dvals, 1e-15)))
}

calc_kl_divergence <- function(N, dgp_fun, model_predict, true_logdens_fun) {
  samp <- dgp_fun(N)
  log_ftrue <- true_logdens_fun(samp)
  fhat      <- model_predict(samp)
  log_fhat  <- log(pmax(fhat, 1e-15))
  mean(log_ftrue - log_fhat)
}

analyze_models <- function(dtrain, dtest = NULL,
                           config, distLibrary,
                           showTest = TRUE) {
  # Nur Demo für d=3-Fall:
  tf <- fit_transformation_forest(dtrain)
  kd <- fit_kernel_density_3d(dtrain)
  pm <- fit_simple_param(dtrain)
  
  tf_dtrain <- tf$predict_joint(dtrain)
  kd_dtrain <- kd$predict_kde(dtrain)
  pm_dtrain <- pm$predict_joint(dtrain)
  
  tf_ll <- logLik_on_data(tf_dtrain)
  kd_ll <- logLik_on_data(kd_dtrain)
  pm_ll <- logLik_on_data(pm_dtrain)
  
  cat("\n[Train] TRFT LogLik =", tf_ll,
      "\n[Train] KDE  LogLik =", kd_ll,
      "\n[Train] PM   LogLik =", pm_ll, "\n")
  
  if(!is.null(dtest) && showTest) {
    tf_dtest <- tf$predict_joint(dtest)
    kd_dtest <- kd$predict_kde(dtest)
    pm_dtest <- pm$predict_joint(dtest)
    cat("\n[Test ] TRFT LogLik =", logLik_on_data(tf_dtest),
        "\n[Test ] KDE  LogLik =", logLik_on_data(kd_dtest),
        "\n[Test ] PM   LogLik =", logLik_on_data(pm_dtest), "\n")
  }
  
  true_logdens_fun <- function(df) {
    compute_logdensity(df, d=3, config=config, distLibrary=distLibrary)
  }
  
  kl_tf <- calc_kl_divergence(
    N=500,
    dgp_fun = function(nn) generate_data(nn, d=3, config=config, distLibrary=distLibrary),
    model_predict = tf$predict_joint,
    true_logdens_fun = true_logdens_fun
  )
  kl_kd <- calc_kl_divergence(
    N=500,
    dgp_fun = function(nn) generate_data(nn, d=3, config=config, distLibrary=distLibrary),
    model_predict = kd$predict_kde,
    true_logdens_fun = true_logdens_fun
  )
  kl_pm <- calc_kl_divergence(
    N=500,
    dgp_fun = function(nn) generate_data(nn, d=3, config=config, distLibrary=distLibrary),
    model_predict = pm$predict_joint,
    true_logdens_fun = true_logdens_fun
  )
  
  cat("\nApprox KL Divergence:",
      "\n   TF  =", kl_tf,
      "\n   KDE =", kl_kd,
      "\n   PM  =", kl_pm, "\n")
  
  invisible(list(
    KL      = c(TF=kl_tf, KDE=kl_kd, PM=kl_pm),
    LLtrain = c(TF=tf_ll, KDE=kd_ll, PM=pm_ll)
  ))
}





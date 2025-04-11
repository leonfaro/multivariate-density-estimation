# 04_analyze.R

analyze_models <- function(dtrain, dtest=NULL) {
  tf <- fit_transformation_forest(dtrain)
  kd <- fit_kernel_density_3d(dtrain)
  pm <- fit_simple_param(dtrain)
  
  tf_dtrain <- tf$predict_joint(dtrain)
  kd_dtrain <- kd$predict_kde(dtrain)
  pm_dtrain <- pm$predict_joint(dtrain)
  
  tf_ll <- logLik_on_data(tf_dtrain)
  kd_ll <- logLik_on_data(kd_dtrain)
  pm_ll <- logLik_on_data(pm_dtrain)
  
  cat("\n[Train] TRFT  LogLik =", tf_ll,
      "\n[Train] KDE LogLik=", kd_ll,
      "\n[Train] PM  LogLik =", pm_ll, "\n")
  
  if(!is.null(dtest)) {
    tf_dtest <- tf$predict_joint(dtest)
    kd_dtest <- kd$predict_kde(dtest)
    pm_dtest <- pm$predict_joint(dtest)
    cat("\n[Test ] TRFT  LogLik =", logLik_on_data(tf_dtest),
        "\n[Test ] KDE LogLik=", logLik_on_data(kd_dtest),
        "\n[Test ] PM  LogLik =", logLik_on_data(pm_dtest), "\n")
  }
  
  kl_tf <- calc_kl_divergence(50, function(n) my_dgp_d(N=n, d=3, models=models, cond=cond), tf$predict_joint)
  kl_kd <- calc_kl_divergence(50, function(n) my_dgp_d(N=n, d=3, models=models, cond=cond), function(dd) kd$predict_kde(dd))
  kl_pm <- calc_kl_divergence(50, function(n) my_dgp_d(N=n, d=3, models=models, cond=cond), pm$predict_joint)
  
  cat("\nApprox KL Divergence:\n",
      "TF =", kl_tf, "\n",
      "KDE=", kl_kd, "\n",
      "PM =", kl_pm, "\n")
  
  invisible(list(
    KL = c(TF = kl_tf, KDE = kl_kd, PM = kl_pm),
    LLtrain = c(TF = tf_ll, KDE = kd_ll, PM = pm_ll)
  ))
}




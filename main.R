source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ks_model.R")
source("04_evaluation.R")
source("EDA.R")

n <- 50
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "beta", parm = function(d) list(shape1 = softplus(d$X2), shape2 = softplus(d$X1))),
  list(distr = "gamma", parm = function(d) list(shape = softplus(d$X3), scale = softplus(d$X2)))
)
perm <- c(3, 4, 1, 2)

#' @export
main <- function() {
  t0 <- proc.time()

  G <- list(
    n = n,
    config = config,
    seed = 42,
    split_ratio = 0.5
  )

  X <- gen_samples(G)                             # Script 2
  S <- train_test_split(X, G$split_ratio, G$seed) # Script 3

  ## Modelle in Originalreihenfolge inklusive Laufzeitmessung
  t_true <- system.time({
    M_TRUE <- fit_TRUE(S$X_tr, S$X_te, G$config)    # Script 5
    baseline_ll <- logL_TRUE_dim(M_TRUE, S$X_te)
  })[["elapsed"]]

  t_trtf <- system.time({
    M_TRTF <- fit_TRTF(S$X_tr, S$X_te, G$config)    # Script models/trtf_model.R
    trtf_ll <- logL_TRTF_dim(M_TRTF, S$X_te)
  })[["elapsed"]]

  t_ks <- system.time({
    M_KS <- fit_KS(S$X_tr, S$X_te, G$config)        # Script models/ks_model.R
    ks_ll <- logL_KS_dim(M_KS, S$X_te)
  })[["elapsed"]]


  tab_normal <- data.frame(
    dim = as.character(seq_along(G$config)),
    distribution = sapply(G$config, `[[`, "distr"),
    logL_baseline = baseline_ll,
    logL_trtf = trtf_ll,
    logL_ks = ks_ll,
    stringsAsFactors = FALSE
  )

  tab_normal <- add_sum_row(tab_normal)

  ## Modelle mit Permutation
  X_tr_p <- S$X_tr[, perm, drop = FALSE]
  X_te_p <- S$X_te[, perm, drop = FALSE]
  colnames(X_tr_p) <- paste0("X", seq_len(ncol(X_tr_p)))
  colnames(X_te_p) <- paste0("X", seq_len(ncol(X_te_p)))

  t_true_p <- system.time({
    M_TRUE_p <- fit_TRUE(X_tr_p, X_te_p, G$config)
    baseline_ll_p <- logL_TRUE_dim(M_TRUE_p, X_te_p)
  })[["elapsed"]]

  t_trtf_p <- system.time({
    M_TRTF_p <- fit_TRTF(X_tr_p, X_te_p, G$config)
    trtf_ll_p <- logL_TRTF_dim(M_TRTF_p, X_te_p)
  })[["elapsed"]]

  t_ks_p <- system.time({
    M_KS_p <- fit_KS(X_tr_p, X_te_p, G$config)
    ks_ll_p <- logL_KS_dim(M_KS_p, X_te_p)
  })[["elapsed"]]


  tab_perm <- data.frame(
    dim = as.character(seq_along(G$config)),
    distribution = sapply(G$config, `[[`, "distr"),
    logL_baseline = baseline_ll_p,
    logL_trtf = trtf_ll_p,
    logL_ks = ks_ll_p,
    stringsAsFactors = FALSE
  )

  tab_perm <- add_sum_row(tab_perm)

  ## Log-Dichten fuer Plot (Originalreihenfolge)
  ld_base <- rowSums(vapply(seq_along(G$config), function(k) {
    .log_density_vec(S$X_te[, k], G$config[[k]]$distr, M_TRUE$theta[[k]])
  }, numeric(nrow(S$X_te))))
  ld_trtf <- as.numeric(predict(M_TRTF, S$X_te, type = "logdensity"))
  ld_ks   <- as.numeric(predict(M_KS, S$X_te, type = "logdensity"))

  runtime <- round((proc.time() - t0)[["elapsed"]], 2)
  message(paste0("Laufzeit: ", runtime, " Sekunden"))




  ## Log-Dichten fuer Plot (Permutation)
  ld_base_p <- rowSums(vapply(seq_along(G$config), function(k) {
    .log_density_vec(X_te_p[, k], G$config[[k]]$distr, M_TRUE_p$theta[[k]])
  }, numeric(nrow(X_te_p))))
  ld_trtf_p <- as.numeric(predict(M_TRTF_p, X_te_p, type = "logdensity"))
  ld_ks_p   <- as.numeric(predict(M_KS_p,   X_te_p, type = "logdensity"))

  scatter_data <- list(
    ld_base   = ld_base,
    ld_trtf   = ld_trtf,
    ld_ks     = ld_ks,
    ld_base_p = ld_base_p,
    ld_trtf_p = ld_trtf_p,
    ld_ks_p   = ld_ks_p
  )


  vec_normal <- c(true = t_true, trtf = t_trtf, ks = t_ks)
  vec_perm   <- c(true = t_true_p, trtf = t_trtf_p, ks = t_ks_p)
  kbl_tab <- combine_logL_tables(tab_normal, tab_perm,
                                 vec_normal, vec_perm)

  create_EDA_report(S$X_tr, G$config,
                    scatter_data = scatter_data,
                    table_kbl = kbl_tab)

  message(paste0("N = ", G$n))
  message("Beste Hyperparameter fuer TRTF:")
  print(M_TRTF$best_cfg)
  print(attr(kbl_tab, "tab_data"))
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    make_plot <- function(x, y, ttl) {
      ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
        geom_point(color = "steelblue", size = 1) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        labs(title = ttl, x = "True log-Dichte", y = "Geschaetzte log-Dichte") +
        theme_minimal()
    }
    plots <- with(scatter_data, list(
      make_plot(ld_base,   ld_trtf,   "TRTF"),
      make_plot(ld_base_p, ld_trtf_p, "TRTF (Perm.)"),
      make_plot(ld_base,   ld_ks,     "KS"),
      make_plot(ld_base_p, ld_ks_p,   "KS (Perm.)")
    ))
    print(gridExtra::grid.arrange(grobs = plots, ncol = 2))
  }

  res <- kbl_tab
  attr(res, "tab_data") <- attr(kbl_tab, "tab_data")
  res
}



if (sys.nframe() == 0L) {
  main()
}


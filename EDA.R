#' Exploratory Data Analysis for Generated Samples
#'
#' Creates tables and Scatter-Plots der geschaetzten gegen die wahren
#' Log-Dichten.  Gibt Plots und Tabellen als Liste zurueck.
#'
#' @param X matrix of generated samples
#' @param cfg configuration list as in `gen_samples`
#' @param output_file bislang Pfad zu einer PDF-Datei (wird ignoriert)
#' @return Liste mit Elementen `plots`, `param_plots` und `table`

create_param_plots <- function(param_list) {
  plots <- list()
  idx <- 1L
  for (k in seq_along(param_list)) {
    df <- param_list[[k]]
    if (k == 1 || is.null(df)) next
    for (nm in names(df)) {
      plt <- ggplot(df, aes_string(x = nm)) +
        geom_histogram(bins = 30, fill = "steelblue", color = "black") +
        labs(title = paste0("X", k, ": ", nm), x = nm, y = "Haeufigkeit") +
        theme_minimal()
      plots[[idx]] <- plt
      idx <- idx + 1L
    }
  }
  plots
}

create_EDA_report <- function(X, cfg, output_file = "eda_report.pdf",
                              scatter_data = NULL,
                              table_kbl = NULL,
                              param_list = NULL) {
  if (!requireNamespace("gridExtra", quietly = TRUE))
    install.packages("gridExtra", repos = "https://cloud.r-project.org")



  make_plot <- function(x, y, ttl) {
    ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
      geom_point(color = "steelblue", size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      labs(title = ttl, x = "True log-Dichte", y = "Geschaetzte log-Dichte") +
      theme_minimal()
  }

  plots <- NULL
  param_plots <- NULL
  if (!is.null(scatter_data)) {
    with(scatter_data, {
      plots <<- list(
        make_plot(ld_base,   ld_trtf,   "TRTF"),
        make_plot(ld_base_p, ld_trtf_p, "TRTF (Perm.)"),
        make_plot(ld_base,   ld_ks,     "KS"),
        make_plot(ld_base_p, ld_ks_p,   "KS (Perm.)")
      )
    })
  }
  if (!is.null(param_list))
    param_plots <- create_param_plots(param_list)

  if (!is.null(table_kbl)) {
    tab_data <- attr(table_kbl, "tab_data")
    num_cols <- vapply(tab_data, is.numeric, logical(1))
    for (col in names(tab_data)[num_cols]) {
      last_val <- as.numeric(tab_data[nrow(tab_data), col])
      tab_data[-nrow(tab_data), col] <- sprintf("%.3f", as.numeric(tab_data[-nrow(tab_data), col]))
      tab_data[nrow(tab_data), col]  <- sprintf("%.0f", last_val)
    }
  }
  if (!is.null(plots))
    gridExtra::grid.arrange(grobs = plots, ncol = 2)
  if (!is.null(param_plots)) {
    for (pg in param_plots)
      gridExtra::grid.arrange(pg)
  }
  list(plots = plots, param_plots = param_plots, table = table_kbl)
}

#' VollstÃ¤ndigen Analyseablauf ausfÃ¼hren
#'
#' Diese Funktion Ã¼bernimmt die Hauptlogik aus `main()` und erzeugt sowohl
#' die Tabellen als auch die Plots der explorativen Datenanalyse. Alle
#' benÃ¶tigten Einstellungen werden als Argumente Ã¼bergeben.
#'
#' @param n      StichprobengrÃ¶Ãe
#' @param config Listenkonfiguration der Zielverteilungen
#' @param perm   Permutationsvektor der Spalten
#' @return `knitr_kable`-Objekt mit zusammengefassten Log-Likelihoods
#' @export
run_pipeline <- function(n, config, perm) {
  t0 <- proc.time()

  G <- list(
    n = n,
    config = config,
    seed = 42,
    split_ratio = 0.5
  )

  X_dat <- gen_samples(G, return_params = TRUE)
  X <- X_dat$X
  param_list <- X_dat$params
  S <- train_test_split(X, G$split_ratio, G$seed)

  M_TRUE <- fit_TRUE(S$X_tr, S$X_te, G$config)
  t_true <- system.time({
    baseline_ll <- logL_TRUE_dim(M_TRUE, S$X_te)
  })[["elapsed"]]

  M_TRTF <- fit_TRTF(S$X_tr, S$X_te, G$config)
  t_trtf <- system.time({
    trtf_ll <- logL_TRTF_dim(M_TRTF, S$X_te)
  })[["elapsed"]]

  M_KS <- fit_KS(S$X_tr, S$X_te, G$config)
  t_ks <- system.time({
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

  X_tr_p <- S$X_tr[, perm, drop = FALSE]
  X_te_p <- S$X_te[, perm, drop = FALSE]
  colnames(X_tr_p) <- paste0("X", seq_len(ncol(X_tr_p)))
  colnames(X_te_p) <- paste0("X", seq_len(ncol(X_te_p)))

  M_TRUE_p <- fit_TRUE(X_tr_p, X_te_p, G$config)
  t_true_p <- system.time({
    baseline_ll_p <- logL_TRUE_dim(M_TRUE_p, X_te_p)
  })[["elapsed"]]

  M_TRTF_p <- fit_TRTF(X_tr_p, X_te_p, G$config)
  t_trtf_p <- system.time({
    trtf_ll_p <- logL_TRTF_dim(M_TRTF_p, X_te_p)
  })[["elapsed"]]

  M_KS_p <- fit_KS(X_tr_p, X_te_p, G$config)
  t_ks_p <- system.time({
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

  ld_base <- rowSums(vapply(seq_along(G$config), function(k) {
    .log_density_vec(S$X_te[, k], G$config[[k]]$distr, M_TRUE$theta[[k]])
  }, numeric(nrow(S$X_te))))
  ld_trtf <- as.numeric(predict(M_TRTF, S$X_te, type = "logdensity"))
  ld_ks   <- as.numeric(predict(M_KS, S$X_te, type = "logdensity"))

  runtime <- round((proc.time() - t0)[["elapsed"]], 2)
  message(paste0("Laufzeit: ", runtime, " Sekunden"))

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
                    table_kbl = kbl_tab,
                    param_list = param_list)

  message(paste0("N = ", G$n))
  message("Beste Hyperparameter fuer TRTF:")
  print(M_TRTF$best_cfg)
  print(attr(kbl_tab, "tab_data"))

  res <- kbl_tab
  attr(res, "tab_data") <- attr(kbl_tab, "tab_data")
  res
}

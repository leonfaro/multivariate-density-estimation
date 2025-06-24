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
      tmp <- tempfile(fileext = ".png")
      png(tmp)
      hist(df[[nm]], breaks = 30, col = "steelblue", border = "black",
           main = paste0("X", k, ": ", nm), xlab = nm)
      plots[[idx]] <- recordPlot()
      dev.off()
      unlink(tmp)
      idx <- idx + 1L
    }
  }
  plots
}

create_EDA_report <- function(X, cfg, output_file = "eda_report.pdf",
                              scatter_data = NULL,
                              table_kbl = NULL,
                              param_list = NULL) {



  make_plot <- function(x, y, ttl) {
    tmp <- tempfile(fileext = ".png")
    png(tmp)
    plot(x, y, main = ttl, xlab = "True log-Dichte",
         ylab = "Geschaetzte log-Dichte",
         col = "steelblue", pch = 16)
    abline(a = 0, b = 1, lty = 2)
    plt <- recordPlot()
    dev.off()
    unlink(tmp)
    plt
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
    for (pg in plots) replayPlot(pg)
  if (!is.null(param_plots))
    for (pg in param_plots) replayPlot(pg)
  list(plots = plots, param_plots = param_plots, table = table_kbl)
}

#' Daten vorbereiten und in Trainings-/Testsets aufteilen
#'
#' @param n Stichprobengroesse
#' @param config Listenkonfiguration
#' @param perm Permutationsvektor
#' @return Liste mit Elementen `G`, `S`, `S_perm` und `param_list`
#' @export
prepare_data <- function(n, config, perm) {
  G <- list(n = n, config = config, seed = 42, split_ratio = 0.5)
  X_dat <- gen_samples(G, return_params = TRUE)
  S <- train_test_split(X_dat$X, G$split_ratio, G$seed)

  X_tr_p <- S$X_tr[, perm, drop = FALSE]
  X_te_p <- S$X_te[, perm, drop = FALSE]
  colnames(X_tr_p) <- paste0("X", seq_len(ncol(X_tr_p)))
  colnames(X_te_p) <- paste0("X", seq_len(ncol(X_te_p)))
  S_perm <- list(X_tr = X_tr_p, X_te = X_te_p)

  list(G = G, S = S, S_perm = S_perm, param_list = X_dat$params)
}

#' TRUE, TRTF und KS Modelle fitten
#'
#' @param S Liste mit `X_tr` und `X_te`
#' @param config Konfiguration der Zielverteilungen
#' @return Liste `models`, `ll` und `times`
#' @export
fit_models <- function(S, config) {
  M_TRUE <- fit_TRUE(S$X_tr, S$X_te, config)
  t_true <- system.time({
    ll_true <- logL_TRUE_dim(M_TRUE, S$X_te)
  })[["elapsed"]]

  M_TRTF <- fit_TRTF(S$X_tr, S$X_te, config)
  t_trtf <- system.time({
    ll_trtf <- logL_TRTF_dim(M_TRTF, S$X_te)
  })[["elapsed"]]

  M_KS <- fit_KS(S$X_tr, S$X_te, config)
  t_ks <- system.time({
    ll_ks <- logL_KS_dim(M_KS, S$X_te)
  })[["elapsed"]]

  list(models = list(true = M_TRUE, trtf = M_TRTF, ks = M_KS),
       ll = list(true = ll_true, trtf = ll_trtf, ks = ll_ks),
       times = c(true = t_true, trtf = t_trtf, ks = t_ks))
}

#' Log-Likelihood-Tabellen erzeugen
#'
#' @param models Rueckgabe von `fit_models`
#' @param config Konfiguration
#' @return Datenrahmen mit Summenzeile
#' @export
calc_loglik_tables <- function(models, config) {
  tab <- data.frame(
    dim = as.character(seq_along(config)),
    distribution = sapply(config, `[[`, "distr"),
    logL_baseline = models$ll$true,
    logL_trtf = models$ll$trtf,
    logL_ks = models$ll$ks,
    stringsAsFactors = FALSE
  )
  add_sum_row(tab)
}

#' Daten fuer Scatter-Plots erzeugen
#'
#' @param models Rueckgabe von `fit_models`
#' @param S Liste mit `X_te`
#' @return Liste mit Logdichte-Vektoren
#' @export
make_scatter_data <- function(models, S) {
  G_cfg <- models$models$true$config
  ld_base <- rowSums(vapply(seq_along(G_cfg), function(k) {
    .log_density_vec(S$X_te[, k], G_cfg[[k]]$distr,
                     models$models$true$theta[[k]])
  }, numeric(nrow(S$X_te))))
  ld_trtf <- as.numeric(predict(models$models$trtf, S$X_te,
                                type = "logdensity"))
  ld_ks <- as.numeric(predict(models$models$ks, S$X_te,
                              type = "logdensity"))
  list(ld_base = ld_base, ld_trtf = ld_trtf, ld_ks = ld_ks)
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
#' @return Liste mit Daten, Modellen, Tabellen und Streudaten
#' @export
run_pipeline <- function(n, config, perm) {
  t0 <- proc.time()
  prep <- prepare_data(n, config, perm)

  mods_norm <- fit_models(prep$S, config)
  tab_normal <- calc_loglik_tables(mods_norm, config)
  scat_norm <- make_scatter_data(mods_norm, prep$S)

  mods_perm <- fit_models(prep$S_perm, config)
  tab_perm <- calc_loglik_tables(mods_perm, config)
  scat_perm <- make_scatter_data(mods_perm, prep$S_perm)

  scatter_data <- c(scat_norm,
                    setNames(scat_perm, paste0(names(scat_perm), "_p")))

  kbl_tab <- combine_logL_tables(tab_normal, tab_perm,
                                 mods_norm$times, mods_perm$times)

  create_EDA_report(prep$S$X_tr, config,
                    scatter_data = scatter_data,
                    table_kbl = kbl_tab,
                    param_list = prep$param_list)

  runtime <- round((proc.time() - t0)[["elapsed"]], 2)
  message(paste0("Laufzeit: ", runtime, " Sekunden"))

  res <- list(
    data = prep,
    models = list(normal = mods_norm$models, perm = mods_perm$models),
    tables = list(normal = tab_normal, perm = tab_perm, combined = kbl_tab),
    scatter_data = scatter_data,
    times = list(normal = mods_norm$times, perm = mods_perm$times)
  )
  res
}

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

#' Daten vorbereiten
#'
#' Generiert Stichproben gemäß `config`, splittet sie in Trainings-
#' und Testdaten und legt eine permutierte Variante an.
#'
#' @param n Stichprobengröße
#' @param config Konfigurationsliste
#' @param perm Permutationsvektor
#' @return Liste mit Trainings- und Testdaten sowie Parameterhistorie
prepare_data <- function(n, config, perm) {
  G <- list(n = n, config = config, seed = 42, split_ratio = 0.5)
  X_dat <- gen_samples(G, return_params = TRUE)
  S <- train_test_split(X_dat$X, G$split_ratio, G$seed)
  X_tr_p <- S$X_tr[, perm, drop = FALSE]
  X_te_p <- S$X_te[, perm, drop = FALSE]
  colnames(X_tr_p) <- paste0("X", seq_len(ncol(X_tr_p)))
  colnames(X_te_p) <- paste0("X", seq_len(ncol(X_te_p)))
  list(X_tr = S$X_tr, X_te = S$X_te,
       X_tr_p = X_tr_p, X_te_p = X_te_p,
       param_list = X_dat$params, config = config)
}

#' Modelle fitten und Laufzeiten messen
#'
#' Passt TRUE, TRTF und KS an normale und permutierte Daten an.
#'
#' @param S Liste aus `prepare_data`
#' @param config Konfigurationsliste
#' @return Liste mit Modellen und Laufzeiten
fit_models <- function(S, config) {
  M_TRUE <- fit_TRUE(S$X_tr, S$X_te, config)
  t_true <- system.time({ logL_TRUE_dim(M_TRUE, S$X_te) })[["elapsed"]]
  M_TRTF <- fit_TRTF(S$X_tr, S$X_te, config)
  t_trtf <- system.time({ logL_TRTF_dim(M_TRTF, S$X_te) })[["elapsed"]]
  M_KS <- fit_KS(S$X_tr, S$X_te, config)
  t_ks <- system.time({ logL_KS_dim(M_KS, S$X_te) })[["elapsed"]]

  M_TRUE_p <- fit_TRUE(S$X_tr_p, S$X_te_p, config)
  t_true_p <- system.time({ logL_TRUE_dim(M_TRUE_p, S$X_te_p) })[["elapsed"]]
  M_TRTF_p <- fit_TRTF(S$X_tr_p, S$X_te_p, config)
  t_trtf_p <- system.time({ logL_TRTF_dim(M_TRTF_p, S$X_te_p) })[["elapsed"]]
  M_KS_p <- fit_KS(S$X_tr_p, S$X_te_p, config)
  t_ks_p <- system.time({ logL_KS_dim(M_KS_p, S$X_te_p) })[["elapsed"]]

  list(
    normal = list(true = M_TRUE, trtf = M_TRTF, ks = M_KS),
    perm   = list(true = M_TRUE_p, trtf = M_TRTF_p, ks = M_KS_p),
    runtime_normal = c(true = t_true, trtf = t_trtf, ks = t_ks),
    runtime_perm   = c(true = t_true_p, trtf = t_trtf_p, ks = t_ks_p)
  )
}

#' Log-Likelihood-Tabellen berechnen
#'
#' Erstellt Tabellen für originale und permutierte Reihenfolge.
#'
#' @param models Ergebnis von `fit_models`
#' @param S Liste aus `prepare_data`
#' @return Liste mit Tabellen `tab_normal` und `tab_perm`
calc_loglik_tables <- function(models, S) {
  cfg <- S$config
  baseline_ll <- logL_TRUE_dim(models$normal$true, S$X_te)
  trtf_ll <- logL_TRTF_dim(models$normal$trtf, S$X_te)
  ks_ll <- logL_KS_dim(models$normal$ks, S$X_te)

  tab_normal <- data.frame(
    dim = as.character(seq_along(cfg)),
    distribution = sapply(cfg, `[[`, "distr"),
    logL_baseline = baseline_ll,
    logL_trtf = trtf_ll,
    logL_ks = ks_ll,
    stringsAsFactors = FALSE
  )
  tab_normal <- add_sum_row(tab_normal)

  baseline_ll_p <- logL_TRUE_dim(models$perm$true, S$X_te_p)
  trtf_ll_p <- logL_TRTF_dim(models$perm$trtf, S$X_te_p)
  ks_ll_p <- logL_KS_dim(models$perm$ks, S$X_te_p)

  tab_perm <- data.frame(
    dim = as.character(seq_along(cfg)),
    distribution = sapply(cfg, `[[`, "distr"),
    logL_baseline = baseline_ll_p,
    logL_trtf = trtf_ll_p,
    logL_ks = ks_ll_p,
    stringsAsFactors = FALSE
  )
  tab_perm <- add_sum_row(tab_perm)

  list(tab_normal = tab_normal, tab_perm = tab_perm)
}

#' Daten für Scatterplots zusammenstellen
#'
#' Gibt die wahren und geschätzten Log-Dichten zurück.
#'
#' @param models Ergebnis von `fit_models`
#' @param S Liste aus `prepare_data`
#' @return Liste mit logDichte-Vektoren
make_scatter_data <- function(models, S) {
  cfg <- S$config
  ld_base <- rowSums(vapply(seq_along(cfg), function(k) {
    .log_density_vec(S$X_te[, k], cfg[[k]]$distr,
                     models$normal$true$theta[[k]])
  }, numeric(nrow(S$X_te))))
  ld_trtf <- as.numeric(predict(models$normal$trtf, S$X_te,
                                type = "logdensity"))
  ld_ks <- as.numeric(predict(models$normal$ks, S$X_te,
                              type = "logdensity"))

  ld_base_p <- rowSums(vapply(seq_along(cfg), function(k) {
    .log_density_vec(S$X_te_p[, k], cfg[[k]]$distr,
                     models$perm$true$theta[[k]])
  }, numeric(nrow(S$X_te_p))))
  ld_trtf_p <- as.numeric(predict(models$perm$trtf, S$X_te_p,
                                  type = "logdensity"))
  ld_ks_p <- as.numeric(predict(models$perm$ks, S$X_te_p,
                                type = "logdensity"))

  list(ld_base = ld_base, ld_trtf = ld_trtf, ld_ks = ld_ks,
       ld_base_p = ld_base_p, ld_trtf_p = ld_trtf_p, ld_ks_p = ld_ks_p)
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
#' @return Liste der Einzelergebnisse
#' @export
run_pipeline <- function(n, config, perm) {
  dat <- prepare_data(n, config, perm)
  models <- fit_models(dat, config)
  tables <- calc_loglik_tables(models, dat)
  scatter <- make_scatter_data(models, dat)
  kbl_tab <- combine_logL_tables(tables$tab_normal, tables$tab_perm,
                                 models$runtime_normal, models$runtime_perm)

  create_EDA_report(dat$X_tr, config,
                    scatter_data = scatter,
                    table_kbl = kbl_tab,
                    param_list = dat$param_list)

  list(data = dat, models = models,
       tables = tables, scatter = scatter,
       kbl_tab = kbl_tab)
}

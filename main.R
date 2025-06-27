source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ks_model.R")
source("04_evaluation.R")

n <- 50
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "beta", parm = function(d) list(shape1 = softplus(d$X2), shape2 = softplus(d$X1))),
  list(distr = "gamma", parm = function(d) list(shape = softplus(d$X3), scale = softplus(d$X2)))
)

run_normal <- function(prep, cfg) {
  mods <- fit_models(prep$S, cfg)
  tab_mean <- calc_loglik_tables(mods, cfg)
  tab_sd   <- calc_loglik_sds(mods, prep$S, cfg)
  tab_fmt  <- format_loglik_table(tab_mean, tab_sd)
  print(tab_fmt)
  hyp <- paste(names(mods$models$trtf$best_cfg),
               mods$models$trtf$best_cfg,
               sep = "=", collapse = ", ")
  cat(sprintf("TRTF Predict normal: %.2f s [%s]\n",
              mods$times["trtf"], hyp))
  list(tab = tab_mean, mods = mods)
}


#' Kombinierte Streuplots anzeigen
#'
#' Vier Streudiagramme werden im Layout 2x2 ausgegeben. Die oberen
#' beiden Plots zeigen die Daten in normaler Reihenfolge, unten folgt
#' die Permutation. Links steht jeweils der TRTF-Vergleich, rechts der
#' Vergleich mit KS.
#'
#' @param scatter_data Liste mit Einträgen `ld_base`, `ld_trtf`,
#'   `ld_ks`, `ld_base_p`, `ld_trtf_p`, `ld_ks_p`
#' @return Ein aufgezeichneter Plot
plot_scatter_matrix <- function(scatter_data) {
  op <- par(mfrow = c(2, 2))
  on.exit(par(op))
  with(scatter_data, {
    plot(ld_trtf, ld_base, xlab = "predicted log-density",
         ylab = "true log-density")
    abline(0, 1)
    plot(ld_ks, ld_base, xlab = "predicted log-density",
         ylab = "true log-density")
    abline(0, 1)
    plot(ld_trtf_p, ld_base_p, xlab = "predicted log-density",
         ylab = "true log-density")
    abline(0, 1)
    plot(ld_ks_p, ld_base_p, xlab = "predicted log-density",
         ylab = "true log-density")
    abline(0, 1)
  })
  recordPlot()
}

#' Starte die komplette Analyse
#'
#' Es werden zwei Log-Likelihood-Tabellen ausgegeben: einmal für die
#' normale Reihenfolge der Dimensionen und einmal für eine vorgegebene
#' Permutation.
#'
#' @export
main <- function() {
  prep <- prepare_data(n, config, seq_along(config))

  res_norm <- run_normal(prep, config)

  invisible(list(normal = res_norm$tab))
}

if (sys.nframe() == 0L) {
  main()
}

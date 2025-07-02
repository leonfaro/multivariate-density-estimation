# Utilities for model evaluation

#' Summen-Zeile an Tabelle anh\u00e4ngen
#'
#' @param tab data frame mit numerischen Spalten
#' @param label Bezeichner f\u00fcr die Summenzeile
#' @return data frame mit zus\u00e4tzlicher Zeile
#' @export
add_sum_row <- function(tab, label = "k") {
  stopifnot(is.data.frame(tab))
  sum_row <- setNames(vector("list", ncol(tab)), names(tab))
  for (nm in names(tab)) {
    if (nm == "dim") {
      sum_row[[nm]] <- label
    } else if (is.numeric(tab[[nm]])) {
      sum_row[[nm]] <- sum(tab[[nm]])
    } else {
      sum_row[[nm]] <- NA
    }
  }
  rbind(tab, as.data.frame(sum_row, stringsAsFactors = FALSE))
}

#' Daten erzeugen und aufteilen
#'
#' @param n Stichprobengr\u00f6\u00dfe
#' @param config Verteilungskonfiguration
#' @param seed Zufallsstartwert
#' @return Liste mit Originaldaten `X` und Splits `S`
#' @export
prepare_data <- function(n, config, seed = 42) {
  X <- Generate_iid_from_config(n, config)
  S <- split_data(X, seed)
  list(X = X, S = S)
}

#' TRUE, TRTF und KS Modelle fitten
#'
#' @param S Liste mit `X_tr`, `X_val`, `X_te`
#' @param config Konfiguration der Zielverteilungen
#' @return Liste `models`, `ll` und `times`
#' @export
fit_models <- function(S, config) {
  M_TRUE <- fit_TRUE(S, config)
  t_true <- system.time({
    ll_true <- logL_TRUE_dim(M_TRUE, S$X_te)
  })[["elapsed"]]

  M_TRTF <- fit_TRTF(S, config)
  t_trtf <- system.time({
    ll_trtf <- logL_TRTF_dim(M_TRTF, S$X_te)
  })[["elapsed"]]

  M_KS <- fit_KS(S, config)
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

#' Standardabweichungen der Log-Likelihoods
#'
#' @param models Rueckgabe von `fit_models`
#' @param S Liste mit `X_te`
#' @param config Konfiguration
#' @return Datenrahmen analog zu `calc_loglik_tables`
#' @export
calc_loglik_sds <- function(models, S, config) {
  K <- length(config)
  sd_true <- numeric(K)
  for (k in seq_len(K)) {
    ll_vec <- .log_density_vec(S$X_te[, k], config[[k]]$distr,
                               models$models$true$theta[[k]])
    sd_true[k] <- sd(-ll_vec)
  }
  ll_trtf <- -predict(models$models$trtf, S$X_te,
                      type = "logdensity_by_dim")
  ll_ks   <- -predict(models$models$ks,  S$X_te,
                      type = "logdensity_by_dim")
  sd_trtf <- apply(ll_trtf, 2, sd)
  sd_ks   <- apply(ll_ks,   2, sd)

  tab <- data.frame(
    dim = as.character(seq_len(K)),
    distribution = sapply(config, `[[`, "distr"),
    logL_baseline = sd_true,
    logL_trtf = sd_trtf,
    logL_ks = sd_ks,
    stringsAsFactors = FALSE
  )
  sum_row <- list(
    dim = "k",
    distribution = "SUM",
    logL_baseline = sqrt(sum(sd_true^2)),
    logL_trtf = sqrt(sum(sd_trtf^2)),
    logL_ks = sqrt(sum(sd_ks^2))
  )
  rbind(tab, as.data.frame(sum_row, stringsAsFactors = FALSE))
}

#' Formatierte Log-Likelihood-Tabelle
#'
#' @param tab_mean Mittelwerte
#' @param tab_sd Standardabweichungen
#' @return Datenrahmen mit Zeichenketten
#' @export
format_loglik_table <- function(tab_mean, tab_sd) {
  tab <- tab_mean
  cols <- c("logL_baseline", "logL_trtf", "logL_ks")
  for (col in cols) {
    tab[[col]] <- sprintf("%.2f ± %.2f",
                          round(tab_mean[[col]], 2),
                          round(2 * tab_sd[[col]], 2))
  }
  tab
}

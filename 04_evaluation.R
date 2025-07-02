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

stderr <- function(values) {
  stats::sd(values) / sqrt(length(values))
}

nats_per_dim <- function(loss, d) {
  loss / d
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
calc_loglik_tables <- function(models, config, X_te) {
  K <- length(config)

  ll_true <- matrix(NA_real_, nrow = nrow(X_te), ncol = K)
  for (k in seq_len(K)) {
    ll_vec <- .log_density_vec(X_te[, k], config[[k]]$distr,
                               models$true$theta[[k]])
    ll_true[, k] <- -ll_vec
  }
  ll_trtf <- -predict(models$trtf, X_te, type = "logdensity_by_dim")
  ll_ks   <- -predict(models$ks,  X_te, type = "logdensity_by_dim")

  mean_true <- colMeans(ll_true)
  se_true   <- apply(ll_true, 2, stderr)
  mean_trtf <- colMeans(ll_trtf)
  se_trtf   <- apply(ll_trtf, 2, stderr)
  mean_ks   <- colMeans(ll_ks)
  se_ks     <- apply(ll_ks,   2, stderr)

  fmt <- function(m, se) sprintf("%.2f ± %.2f", round(m, 2), round(2 * se, 2))

  tab <- data.frame(
    dim = as.character(seq_len(K)),
    distribution = sapply(config, `[[`, "distr"),
    true = fmt(mean_true, se_true),
    trtf = fmt(mean_trtf, se_trtf),
    ks   = fmt(mean_ks, se_ks),
    ttm  = rep(NA_character_, K),
    stringsAsFactors = FALSE
  )

  sum_row <- data.frame(
    dim = "k",
    distribution = "SUM",
    true = fmt(sum(mean_true), sqrt(sum(se_true^2))),
    trtf = fmt(sum(mean_trtf), sqrt(sum(se_trtf^2))),
    ks   = fmt(sum(mean_ks),   sqrt(sum(se_ks^2))),
    ttm  = NA_character_,
    stringsAsFactors = FALSE
  )
  rbind(tab, sum_row)
}

#' Standardabweichungen der Log-Likelihoods
#'
#' @param models Rueckgabe von `fit_models`
#' @param S Liste mit `X_te`
#' @param config Konfiguration
#' @return Datenrahmen analog zu `calc_loglik_tables`
#' @export

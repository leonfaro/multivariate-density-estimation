# Utilities for model evaluation
root_path <- getwd()
if (basename(root_path) == "testthat") {
  root_path <- dirname(dirname(root_path))
}
source(file.path(root_path, "models/ttm_marginal.R"))
source(file.path(root_path, "models/ttm_separable.R"))
source(file.path(root_path, "models/ks_model.R"))
source(file.path(root_path, "models/ttm_cross_term.R"))
source(file.path(root_path, "models/true_joint_model.R"))

#' Summen-Zeile an Tabelle anh\u00e4ngen
#'
#' @param tab data frame mit numerischen Spalten
#' @param label Bezeichner f\u00fcr die Summenzeile
#' @details Fehlende Werte in numerischen Spalten werden bei der Summenbildung ignoriert.
#' @return data frame mit zus\u00e4tzlicher Zeile
#' @export
add_sum_row <- function(tab, label = "k") {
  stopifnot(is.data.frame(tab))
  sum_row <- setNames(vector("list", ncol(tab)), names(tab))
  for (nm in names(tab)) {
    if (nm == "dim") {
      sum_row[[nm]] <- label
    } else if (is.numeric(tab[[nm]])) {
      sum_row[[nm]] <- sum(tab[[nm]], na.rm = TRUE)
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

  M_TTM_cross <- trainCrossTermMap(S)
  t_ttm_cross <- system.time({
    ll_ttm_cross <- -predict(M_TTM_cross$S, S$X_te, type = "logdensity_by_dim")
  })[["elapsed"]]

  list(models = list(true = M_TRUE, trtf = M_TRTF, ks = M_KS,
                     ttm_cross = M_TTM_cross),
       ll = list(true = ll_true, trtf = ll_trtf, ks = ll_ks,
                ttm_cross = ll_ttm_cross),
       times = c(true = t_true, trtf = t_trtf, ks = t_ks,
                 ttm_cross = t_ttm_cross))
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
  ll_true_joint <- - true_joint_logdensity_by_dim(config, X_te)
  if (!is.null(models$ttm)) {
    ll_ttm <- -predict(models$ttm$S, X_te, type = "logdensity_by_dim")
    mean_ttm <- colMeans(ll_ttm)
    se_ttm   <- apply(ll_ttm, 2, stderr)
    total_nll_ttm <- rowSums(ll_ttm)
    se_sum_ttm <- stats::sd(total_nll_ttm) / sqrt(length(total_nll_ttm))
  } else {
    mean_ttm <- rep(NA_real_, K)
    se_ttm   <- rep(NA_real_, K)
    se_sum_ttm <- NA_real_
  }
  if (!is.null(models$ttm_sep)) {
    ll_ttm_sep <- -predict(models$ttm_sep$S, X_te, type = "logdensity_by_dim")
    mean_sep <- colMeans(ll_ttm_sep)
    se_sep   <- apply(ll_ttm_sep, 2, stderr)
    total_nll_sep <- rowSums(ll_ttm_sep)
    se_sum_sep <- stats::sd(total_nll_sep) / sqrt(length(total_nll_sep))
  } else {
    mean_sep <- rep(NA_real_, K)
    se_sep   <- rep(NA_real_, K)
    se_sum_sep <- NA_real_
  }
  if (!is.null(models$ttm_cross)) {
    ll_ttm_cross <- -predict(models$ttm_cross$S, X_te,
                             type = "logdensity_by_dim")
    mean_cross <- colMeans(ll_ttm_cross)
    se_cross   <- apply(ll_ttm_cross, 2, stderr)
    total_nll_cross <- rowSums(ll_ttm_cross)
    se_sum_cross <- stats::sd(total_nll_cross) / sqrt(length(total_nll_cross))
  } else {
    mean_cross <- rep(NA_real_, K)
    se_cross   <- rep(NA_real_, K)
    se_sum_cross <- NA_real_
  }

  mean_true <- colMeans(ll_true)
  se_true   <- apply(ll_true, 2, stderr)
  total_nll_true <- rowSums(ll_true)
  se_sum_true <- stats::sd(total_nll_true) / sqrt(length(total_nll_true))
  mean_true_joint <- colMeans(ll_true_joint)
  se_true_joint   <- apply(ll_true_joint, 2, stderr)
  total_nll_true_joint <- rowSums(ll_true_joint)
  se_sum_true_joint <- stats::sd(total_nll_true_joint) /
    sqrt(length(total_nll_true_joint))
  mean_trtf <- colMeans(ll_trtf)
  se_trtf   <- apply(ll_trtf, 2, stderr)
  total_nll_trtf <- rowSums(ll_trtf)
  se_sum_trtf <- stats::sd(total_nll_trtf) / sqrt(length(total_nll_trtf))
  mean_ks   <- colMeans(ll_ks)
  se_ks     <- apply(ll_ks,   2, stderr)
  total_nll_ks <- rowSums(ll_ks)
  se_sum_ks <- stats::sd(total_nll_ks) / sqrt(length(total_nll_ks))
  fmt <- function(m, se) sprintf("%.2f ± %.2f", round(m, 2), round(2 * se, 2))

  tab <- data.frame(
    dim = as.character(seq_len(K)),
    distribution = sapply(config, `[[`, "distr"),
    true = fmt(mean_true, se_true),
    true_joint = fmt(mean_true_joint, se_true_joint),
    trtf = fmt(mean_trtf, se_trtf),
    ks   = fmt(mean_ks, se_ks),
    ttm  = fmt(mean_ttm, se_ttm),
    ttm_sep = fmt(mean_sep, se_sep),
    ttm_cross = fmt(mean_cross, se_cross),
    stringsAsFactors = FALSE
  )

  sum_row <- data.frame(
    dim = "k",
    distribution = "SUM",
    true = fmt(sum(mean_true), se_sum_true),
    true_joint = fmt(sum(mean_true_joint), se_sum_true_joint),
    trtf = fmt(sum(mean_trtf), se_sum_trtf),
    ks   = fmt(sum(mean_ks),   se_sum_ks),
    ttm  = fmt(sum(mean_ttm),  se_sum_ttm),
    ttm_sep = fmt(sum(mean_sep),  se_sum_sep),
    ttm_cross = fmt(sum(mean_cross), se_sum_cross),
    stringsAsFactors = FALSE
  )
  tab <- rbind(tab, sum_row)

  # rename columns for display (do this as the last step before returning)
  nm <- names(tab)
  nm[nm == "true"] <- "True (marginal)"
  nm[nm == "true_joint"] <- "True (Joint)"
  nm[nm == "trtf"] <- "Random Forest"
  nm[nm == "ttm"]  <- "Marginal Map"
  nm[nm == "ttm_sep"] <- "Separable Map"
  nm[nm == "ttm_cross"] <- "Cross-term Map"
  names(tab) <- nm
  message("Ergebnis (NLL in nats; lower is better)")
  tab
}

#' Standardabweichungen der Log-Likelihoods
#'
#' @param models Rueckgabe von `fit_models`
#' @param S Liste mit `X_te`
#' @param config Konfiguration
#' @return Datenrahmen analog zu `calc_loglik_tables`
#' @export

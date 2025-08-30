# Utilities for model evaluation
root_path <- getwd()
if (basename(root_path) == "testthat") {
  root_path <- dirname(dirname(root_path))
}
source(file.path(root_path, "models/ttm_marginal.R"))
source(file.path(root_path, "models/ttm_separable.R"))
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

#' TRUE, TRTF und Cross-Term-Map fitten
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

    M_TTM_cross <- trainCrossTermMap(S)
    t_ttm_cross <- system.time({
      ll_ttm_cross <- -predict(M_TTM_cross$S, S$X_te, type = "logdensity_by_dim")
    })[["elapsed"]]

    list(models = list(true = M_TRUE, trtf = M_TRTF,
                       ttm_cross = M_TTM_cross),
         ll = list(true = ll_true, trtf = ll_trtf,
                   ttm_cross = ll_ttm_cross),
         times = c(true = t_true, trtf = t_trtf,
                   ttm_cross = t_ttm_cross))
}

#' Log-Likelihood-Tabellen erzeugen
#'
#' @param models Rueckgabe von `fit_models`
#' @param config Konfiguration
#' @return Datenrahmen mit Summenzeile
#' @export
calc_loglik_tables <- function(models, config, X_te, config_canonical = NULL, perm = NULL) {
  K <- length(config)

  ll_true <- matrix(NA_real_, nrow = nrow(X_te), ncol = K)
  for (k in seq_len(K)) {
    ll_vec <- .log_density_vec(X_te[, k], config[[k]]$distr,
                               models$true$theta[[k]])
    ll_true[, k] <- -ll_vec
  }
  ll_trtf <- -predict(models$trtf, X_te, type = "logdensity_by_dim")
  ll_true_joint <- tryCatch({
    if (!is.null(config_canonical) && !is.null(perm)) {
      stopifnot(length(perm) == K)
      # Reorder X_te back to the canonical (generative) order, evaluate, then map to current order
      inv_perm <- integer(K); inv_perm[perm] <- seq_len(K)
      X_te_canon <- X_te[, inv_perm, drop = FALSE]
      ll_can <- - true_joint_logdensity_by_dim(config_canonical, X_te_canon)
      stopifnot(ncol(ll_can) == K)
      ll_can[, perm, drop = FALSE]
    } else {
      - true_joint_logdensity_by_dim(config, X_te)
    }
  }, error = function(e) {
    message("[WARN] True (Joint) not computed for this configuration/permutation: ", e$message)
    matrix(NA_real_, nrow = nrow(X_te), ncol = K)
  })

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
  fmt <- function(m, se) sprintf("%.2f ± %.2f", round(m, 2), round(2 * se, 2))

  tab <- data.frame(
    dim = as.character(seq_len(K)),
    distribution = sapply(config, `[[`, "distr"),
    true = fmt(mean_true, se_true),
    true_joint = fmt(mean_true_joint, se_true_joint),
    trtf = fmt(mean_trtf, se_trtf),
    ttm  = fmt(mean_ttm, se_ttm),
    ttm_sep = fmt(mean_sep, se_sep),
    ttm_cross = fmt(mean_cross, se_cross),
    train_test_policy = rep("train_test_only", K),
    stringsAsFactors = FALSE
  )

  sum_row <- data.frame(
    dim = "k",
    distribution = "SUM",
    true = fmt(sum(mean_true), se_sum_true),
    true_joint = fmt(sum(mean_true_joint), se_sum_true_joint),
    trtf = fmt(sum(mean_trtf), se_sum_trtf),
    ttm  = fmt(sum(mean_ttm),  se_sum_ttm),
    ttm_sep = fmt(sum(mean_sep),  se_sum_sep),
    ttm_cross = fmt(sum(mean_cross), se_sum_cross),
    train_test_policy = "train_test_only",
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
  message("Ergebnis (NLL in nats; lower is better) [train/test only]")
  tab
}

#' Standardabweichungen der Log-Likelihoods
#' 
#' @param models Rueckgabe von `fit_models`
#' @param S Liste mit `X_te`
#' @param config Konfiguration
#' @return Datenrahmen analog zu `calc_loglik_tables`
#' @export

#' Evaluate models on two-moons data
#'
#' Fits missing models, computes negative log-likelihoods (nats) and
#' writes a CSV summary.
#'
#' @param mods optional list of fitted models
#' @param S split structure from `make_halfmoon_splits`
#' @param out_csv_path optional output path
#' @return data frame with NLL metrics
#' @export
eval_halfmoon <- function(mods, S, out_csv_path = NULL) {
  dir.create("results", showWarnings = FALSE)
  N <- nrow(S$X_te)
  K <- ncol(S$X_te)
  need <- c("true", "trtf", "ttm", "ttm_sep", "ttm_cross")
  config_moon <- list(list(distr = "norm"), list(distr = "norm"))
  if (missing(mods) || length(mods) == 0 || !all(need %in% names(mods))) {
    seed <- if (!is.null(S$meta$seed)) as.integer(S$meta$seed) else 42L
    set.seed(seed)
    mods <- list(
      true = fit_TRUE(S, config_moon),
      trtf = fit_TRTF(S, config_moon, seed = seed),
      ttm = trainMarginalMap(S, seed = seed)$S,
      ttm_sep = trainSeparableMap(S, seed = seed)$S,
      ttm_cross = trainCrossTermMap(S, seed = seed)$S
    )
  }
  rows <- list()
  for (m in need) {
    mod <- mods[[m]]
    if (m == "true") {
      source("scripts/true_halfmoon_density.R")
      te_true <- true_logdensity(S$X_te, S, Q = 32L)
      LD <- te_true$by_dim
    } else {
      LD <- predict(mod, S$X_te, "logdensity_by_dim")
    }
    stopifnot(is.matrix(LD), all(dim(LD) == c(N, K)), all(is.finite(LD)))
    LDj <- rowSums(LD)
    stopifnot(length(LDj) == N, all(is.finite(LDj)),
              max(abs(rowSums(LD) - LDj)) < 1e-10)
    nllj <- -LDj
    per <- -colMeans(LD)
    se <- stats::sd(nllj) / sqrt(N)
    rows[[length(rows) + 1]] <- c(
      list(model = m, mean_joint_nll = mean(nllj), se_joint = se,
           train_test_policy = "train_test_only"),
      setNames(as.list(per), paste0("per_dim_nll_", 1:K))
    )
  }
  df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  path <- if (is.null(out_csv_path))
    sprintf("results/nll_halfmoon_seed%03d.csv", as.integer(S$meta$seed))
  else out_csv_path
  write.csv(df, path, row.names = FALSE)
  print(df)
  results_table <<- df
  df
}

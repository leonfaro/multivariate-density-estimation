# Utilities for model evaluation. Assumes `source_repo_modules()` has been
# executed so that model components and helpers are available in the search
# path. Standalone sourcing can call `initialize_repo()` first.

# Local stderr helper (avoid base::stderr connection)
stderr <- function(x) stats::sd(x) / sqrt(length(x))

#' Kolmogorov–Smirnov distance to Uniform(0,1)
#'
#' @param u numeric vector in [0,1]
#' @return KS D statistic as numeric scalar
#' @export
ks_distance_u01 <- function(u) {
  u <- as.numeric(u)
  u <- u[is.finite(u)]
  if (!length(u)) return(NA_real_)
  suppressWarnings(as.numeric(stats::ks.test(u, "punif")$statistic))
}

#' Probability integral transform (PIT) for TTM models by dimension
#'
#' @param model TTM model or fit with $S
#' @param X numeric matrix of observations
#' @return N x K matrix of Uniform(0,1) PIT values (Phi(Z_k))
#' @export
pit_ttm_by_dim <- function(model, X) {
  stopifnot(is.matrix(X))
  m <- if (is.list(model) && !is.null(model$S)) model$S else model
  out <- ttm_forward(m, X)
  stats::pnorm(out$Z)
}

#' PIT for TRUE marginal model by dimension
#'
#' @param M_TRUE object from fit_TRUE
#' @param X numeric matrix
#' @return N x K matrix of Uniform(0,1) PIT values
#' @export
pit_true_marginal_by_dim <- function(M_TRUE, X) {
  stopifnot(is.matrix(X))
  cfg <- M_TRUE$config; th <- M_TRUE$theta
  K <- length(cfg); N <- nrow(X)
  U <- matrix(NA_real_, N, K)
  for (k in seq_len(K)) {
    distr <- cfg[[k]]$distr
    par <- th[[k]]
    xk <- X[, k]
    if (distr == "norm") {
      U[, k] <- stats::pnorm(xk, mean = par[1], sd = par[2])
    } else if (distr == "exp") {
      U[, k] <- stats::pexp(pmax(xk, 0), rate = par[1])
    } else if (distr == "beta") {
      U[, k] <- stats::pbeta(pmin(pmax(xk, 0), 1), shape1 = par[1], shape2 = par[2])
    } else if (distr == "gamma") {
      U[, k] <- stats::pgamma(pmax(xk, 0), shape = par[1], scale = par[2])
    } else stop("Unsupported TRUE marginal distribution: ", distr)
  }
  pmin(pmax(U, .Machine$double.eps), 1 - .Machine$double.eps)
}

#' PIT for TRUE joint (triangular) model by dimension via Rosenblatt
#'
#' @param config configuration list used for true_joint_logdensity_by_dim
#' @param X numeric matrix (canonical order matching config)
#' @return N x K matrix of Uniform(0,1) PIT values
#' @export
pit_true_joint_by_dim <- function(config, X) {
  stopifnot(is.matrix(X))
  N <- nrow(X); K <- ncol(X)
  U <- matrix(NA_real_, N, K)
  prev_names <- paste0("X", seq_len(K))
  for (i in seq_len(N)) {
    x_row <- X[i, ]
    for (k in seq_len(K)) {
      prev <- if (k == 1) data.frame() else {
        df <- as.data.frame(as.list(x_row[seq_len(k - 1)]))
        names(df) <- prev_names[seq_len(k - 1)]
        df
      }
      distr_k <- config[[k]]$distr
      args <- if (is.null(config[[k]]$parm)) list() else config[[k]]$parm(prev)
      # sanitize like in true_joint_model
      if (distr_k == "gamma" && all(c("shape1", "shape2") %in% names(args))) {
        args$shape <- args$shape1; args$scale <- args$shape2; args$shape1 <- NULL; args$shape2 <- NULL
      }
      xk <- x_row[k]
      U[i, k] <- switch(distr_k,
        norm  = stats::pnorm(xk, mean = args$mean %||% 0, sd = (args$sd %||% 1)),
        exp   = stats::pexp(pmax(xk, 0), rate = args$rate %||% 1),
        beta  = stats::pbeta(pmin(pmax(xk, 0), 1), shape1 = args$shape1 %||% 1, shape2 = args$shape2 %||% 1),
        gamma = stats::pgamma(pmax(xk, 0), shape = args$shape %||% 1, scale = args$scale %||% 1),
        stop("Unsupported distribution in TRUE joint: ", distr_k)
      )
    }
  }
  pmin(pmax(U, .Machine$double.eps), 1 - .Machine$double.eps)
}

#' PIT for TRTF by dimension (best-effort)
#'
#' Attempts to query CDF from the BoxCox base for k=1 and from
#' traforest conditionals for k>=2. Returns NA if unavailable.
#'
#' @param model fitted mytrtf model
#' @param X numeric matrix
#' @return N x K matrix with PITs or NA when unsupported
#' @export
pit_trtf_by_dim <- function(model, X) {
  stopifnot(is.matrix(X), inherits(model, "mytrtf"))
  K <- length(model$ymod); N <- nrow(X)
  df_new <- as.data.frame(X); names(df_new) <- paste0("X", seq_len(K))
  U <- matrix(NA_real_, N, K)
  # First dimension: unconditional BoxCox
  q1 <- df_new[["X1"]]
  U[, 1] <- tryCatch({
    as.numeric(predict(model$ymod[[1]], newdata = df_new, type = "cdf", q = q1))
  }, error = function(e) NA_real_)
  # Subsequent dimensions via traforest
  for (k in 2:K) {
    fr <- model$forests[[k - 1L]]
    q <- df_new[[paste0("X", k)]]
    U[, k] <- tryCatch({
      as.numeric(predict(fr, newdata = df_new, type = "cdf", q = q))
    }, error = function(e) NA_real_)
  }
  U
}

#' PIT for copula_np by dimension (mixture marginals)
#'
#' Uses per-class KDE marginals and class priors to build
#' unconditional marginal CDFs; joint PIT is not computed here.
#' @param model copula_np model
#' @param X numeric matrix (2 columns)
#' @return N x 2 matrix of PIT values
#' @export
pit_copula_np_by_dim <- function(model, X) {
  stopifnot(inherits(model, "copula_np"), is.matrix(X), ncol(X) == 2L)
  if (!requireNamespace("kde1d", quietly = TRUE)) stop("copula_np PIT requires 'kde1d'")
  N <- nrow(X)
  U <- matrix(NA_real_, N, 2L)
  if (identical(model$mode, "kde_product")) {
    # fallback KDE product: approximate CDF via empirical ECDF on train
    warning("copula_np: PIT for kde_product uses ECDF approximation")
    return(apply(X, 2, function(x) {
      u <- rank(x, ties.method = "average") / (length(x) + 1)
      pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
    }))
  }
  classes <- as.character(model$classes)
  priors <- model$priors
  for (j in seq_along(classes)) {
    yy <- classes[[j]]; comp <- model$by_class[[yy]]
    F1 <- pmin(pmax(kde1d::pkde1d(X[, 1], comp$m1), comp$eps), 1 - comp$eps)
    F2 <- pmin(pmax(kde1d::pkde1d(X[, 2], comp$m2), comp$eps), 1 - comp$eps)
    w <- priors[yy]
    if (!is.finite(w)) w <- 1 / length(classes)
    if (j == 1L) { U[, 1] <- w * F1; U[, 2] <- w * F2 } else { U[, 1] <- U[, 1] + w * F1; U[, 2] <- U[, 2] + w * F2 }
  }
  pmin(pmax(U, .Machine$double.eps), 1 - .Machine$double.eps)
}

#' Compute KS distances (per-dim and median) for available models
#'
#' @param mods list of fitted models as in halfmoon/4d pipelines
#' @param S split list (needs X_te)
#' @param config optional config (for TRUE joint)
#' @return data.frame with columns model, ks_per_dim (list-column), ks_median
#' @export
calc_ks_summary <- function(mods, S, config = NULL) {
  X <- S$X_te
  out <- list()
  add <- function(name, U) {
    if (is.null(U) || !is.matrix(U)) return()
    K <- ncol(U)
    ks <- vapply(seq_len(K), function(k) ks_distance_u01(U[, k]), numeric(1))
    out[[length(out) + 1]] <<- data.frame(model = name,
                                         ks_median = stats::median(ks, na.rm = TRUE),
                                         stringsAsFactors = FALSE)
    attr(out[[length(out)]], "ks_per_dim") <<- ks
  }
  if (!is.null(mods$true)) {
    add("True (marginal)", pit_true_marginal_by_dim(mods$true, X))
  }
  if (!is.null(config)) {
    add("True (Joint)", pit_true_joint_by_dim(config, X))
  }
  if (!is.null(mods$trtf)) {
    add("Random Forest", tryCatch(pit_trtf_by_dim(mods$trtf, X), error = function(e) NULL))
  }
  if (!is.null(mods$ttm)) {
    add("Marginal Map", pit_ttm_by_dim(mods$ttm, X))
  }
  if (!is.null(mods$ttm_sep)) {
    add("Separable Map", pit_ttm_by_dim(mods$ttm_sep, X))
  }
  if (!is.null(mods$copula_np)) {
    if (ncol(X) == 2L) add("Copula NP", tryCatch(pit_copula_np_by_dim(mods$copula_np, X), error = function(e) NULL))
  }
  if (!length(out)) return(data.frame())
  df <- do.call(rbind, out)
  ks_list <- lapply(out, function(d) attr(d, "ks_per_dim"))
  df$ks_per_dim <- ks_list
  df
}

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

    list(models = list(true = M_TRUE, trtf = M_TRTF),
         ll = list(true = ll_true, trtf = ll_trtf),
         times = c(true = t_true, trtf = t_trtf))
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
  ll_trtf <- tryCatch({
    -predict(models$trtf, X_te, type = "logdensity_by_dim")
  }, error = function(e) {
    message("[WARN] TRTF not available or failed: ", e$message)
    matrix(NA_real_, nrow = nrow(X_te), ncol = K)
  })
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
    ll_ttm <- tryCatch({
      -predict_ttm(models$ttm$S, X_te, type = "logdensity_by_dim")
    }, error = function(e) {
      -predict(models$ttm$S, X_te, type = "logdensity_by_dim")
    })
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
    ll_ttm_sep <- tryCatch({
      -predict_ttm(models$ttm_sep$S, X_te, type = "logdensity_by_dim")
    }, error = function(e) {
      -predict(models$ttm_sep$S, X_te, type = "logdensity_by_dim")
    })
    mean_sep <- colMeans(ll_ttm_sep)
    se_sep   <- apply(ll_ttm_sep, 2, stderr)
    total_nll_sep <- rowSums(ll_ttm_sep)
    se_sum_sep <- stats::sd(total_nll_sep) / sqrt(length(total_nll_sep))
  } else {
    mean_sep <- rep(NA_real_, K)
    se_sep   <- rep(NA_real_, K)
    se_sum_sep <- NA_real_
  }
  # Cross-term removed

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
  
  names(tab) <- nm
  message("Ergebnis (NLL; lower is better) [train/test only]")
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
#' Fits missing models, computes negative log-likelihoods, and
#' writes a CSV summary.
#'
#' @param mods optional list of fitted models
#' @param S split structure from `make_halfmoon_splits`
#' @param out_csv_path optional output path
#' @return data frame with NLL metrics
#' @export
eval_halfmoon <- function(mods, S, out_csv_path = NULL) {
  no_artifacts <- isTRUE(getOption("mde.no_artifacts", FALSE)) || identical(Sys.getenv("NO_ARTIFACTS", "0"), "1")
  if (!no_artifacts) {
    dir.create("results", showWarnings = FALSE)
  }
  N <- nrow(S$X_te)
  K <- ncol(S$X_te)
  # Rows to evaluate; include both oracle variants
  need_rows <- c("true_unconditional", "true_conditional", "trtf", "ttm", "ttm_sep", "copula_np")
  # Models that require fitting (the oracle rows are computed directly)
  need_mods <- c("trtf", "ttm", "ttm_sep", "copula_np")
  config_moon <- list(list(distr = "norm"), list(distr = "norm"))
  if (missing(mods) || length(mods) == 0 || !all(need_mods %in% names(mods))) {
    seed <- if (!is.null(S$meta$seed)) as.integer(S$meta$seed) else 42L
    set.seed(seed)
    mods <- list(
      trtf = tryCatch(fit_TRTF(S, config_moon, seed = seed), error = function(e) NULL),
      ttm = trainMarginalMap(S, seed = seed)$S,
      ttm_sep = trainSeparableMap(S, seed = seed)$S,
      copula_np = fit_copula_np(S, seed = seed)
    )
  }
  rows <- list()
  for (m in need_rows) {
    mod <- mods[[m]]
    if (m == "true_unconditional") {
      source(file.path(root_path, "experiments", "halfmoon", "true_density.R"))
      te_true <- true_logdensity(S$X_te, S, Q = 32L)
      LD <- te_true$by_dim
      # Invariants: shape, finiteness, and by-dim sums equal joint
      stopifnot(is.matrix(LD), all(dim(LD) == c(N, K)), all(is.finite(LD)))
      stopifnot(length(te_true$joint) == N, all(is.finite(te_true$joint)))
      stopifnot(max(abs(rowSums(LD) - te_true$joint)) <= 1e-12)
    } else if (m == "true_conditional") {
      source(file.path(root_path, "experiments", "halfmoon", "true_density.R"))
      te_cond <- true_logdensity_conditional(S$X_te, S$y_te, S, Q = 32L)
      LD <- te_cond$by_dim
      # Invariants: shape, finiteness, and by-dim sums equal joint
      stopifnot(is.matrix(LD), all(dim(LD) == c(N, K)), all(is.finite(LD)))
      stopifnot(length(te_cond$joint) == N, all(is.finite(te_cond$joint)))
      stopifnot(max(abs(rowSums(LD) - te_cond$joint)) <= 1e-12)
    } else if (m == "copula_np") {
      LD <- predict(mod, S$X_te, type = "logdensity_by_dim")
    } else {
      if (is.null(mod)) {
        LD <- matrix(NA_real_, N, K)
      } else {
        LD <- predict(mod, S$X_te, "logdensity_by_dim")
      }
    }
    stopifnot(is.matrix(LD), all(dim(LD) == c(N, K)), all(is.finite(LD) | is.na(LD)))
    LDj <- rowSums(LD)
    stopifnot(length(LDj) == N, all(is.finite(LDj) | is.na(LDj)))
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
  if (!no_artifacts) {
    path <- if (is.null(out_csv_path))
      sprintf("results/nll_halfmoon_seed%03d.csv", as.integer(S$meta$seed))
    else out_csv_path
    write.csv(df, path, row.names = FALSE)
  }
  print(df)
  results_table <<- df
  df
}

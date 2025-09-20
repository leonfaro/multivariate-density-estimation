if (!exists("initialize_repo")) {
  source(file.path(getwd(), "R", "loader.R"))
}
root_path <- initialize_repo()

RUN_CONFIG <- list(
  enabled = FALSE,
  train_csv = "/Users/leonkiafaro/Documents/Masterthesis/data_power/power_train_head_1pct.csv",
  test_csv  = "/Users/leonkiafaro/Documents/Masterthesis/data_power/power_test_head_1pct.csv",
  n_train = 250L,
  n_test  = 250L,
  models = c("trtf","ttm","ttm_sep","copula_np")
)

# Small helpers
stderr <- function(x) stats::sd(x) / sqrt(length(x))

`%||%` <- function(a, b) if (is.null(a)) b else a

read_numeric_matrix <- function(path, n_rows = NA_integer_) {
  if (!file.exists(path)) stop(sprintf("CSV not found: %s", path))
  df <- utils::read.csv(path, header = TRUE, check.names = TRUE)
  # Keep only numeric columns present
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  if (!length(num_cols)) stop(sprintf("No numeric columns in %s", basename(path)))
  df <- df[, num_cols, drop = FALSE]
  if (is.finite(n_rows) && !is.na(n_rows) && n_rows > 0L) {
    n_keep <- min(n_rows, nrow(df))
    df <- utils::head(df, n_keep)
  }
  X <- as.matrix(df)
  # Drop rows with NA/Inf
  ok <- stats::complete.cases(X) & apply(is.finite(X), 1, all)
  if (!all(ok)) X <- X[ok, , drop = FALSE]
  storage.mode(X) <- "double"
  X
}

parse_args <- function() {
  out <- list()
  for (a in commandArgs(trailingOnly = TRUE)) {
    if (grepl("^--n_train=", a)) out$n_train <- as.integer(sub("^--n_train=", "", a))
    if (grepl("^--n_test=", a))  out$n_test  <- as.integer(sub("^--n_test=",  "", a))
    if (grepl("^--train_csv=", a)) out$train_csv <- sub("^--train_csv=", "", a)
    if (grepl("^--test_csv=", a))  out$test_csv  <- sub("^--test_csv=",  "", a)
    if (grepl("^--models=", a))    out$models    <- sub("^--models=",    "", a)
    if (grepl("^--seed=", a))      out$seed      <- as.integer(sub("^--seed=", "", a))
  }
  out
}

want_model <- function(name, selected) {
  if (is.null(selected) || !length(selected)) return(TRUE)
  name %in% selected
}

main_power <- function() {
  args <- parse_args()
  if (is.list(RUN_CONFIG) && isTRUE(RUN_CONFIG$enabled)) {
    args$train_csv <- RUN_CONFIG$train_csv
    args$test_csv  <- RUN_CONFIG$test_csv
    args$n_train   <- RUN_CONFIG$n_train
    args$n_test    <- RUN_CONFIG$n_test
    args$models    <- paste(RUN_CONFIG$models, collapse = ",")
  }
  seed <- if (!is.null(args$seed)) args$seed else as.integer(Sys.getenv("SEED", 42))
  set.seed(seed)

  # Default CSV paths provided by user
  train_csv <- args$train_csv %||% Sys.getenv("TRAIN_CSV", "/Users/leonkiafaro/Documents/Masterthesis/data_power/power_train_head_0.1pct.csv")
  test_csv  <- args$test_csv  %||% Sys.getenv("TEST_CSV",  "/Users/leonkiafaro/Documents/Masterthesis/data_power/power_test_head_0.1pct.csv")

  # Number of rows to use (first N rows of each set)
  n_train <- if (!is.null(args$n_train)) args$n_train else suppressWarnings(as.integer(Sys.getenv("N_TRAIN_OVERRIDE", "NA")))
  n_test  <- if (!is.null(args$n_test))  args$n_test  else suppressWarnings(as.integer(Sys.getenv("N_TEST_OVERRIDE",  "NA")))
  if (!is.finite(n_train)) n_train <- NA_integer_
  if (!is.finite(n_test))  n_test  <- NA_integer_

  # Select models to run (comma-separated), default: run all
  models_env <- Sys.getenv("MODELS", "")
  models_sel <- trimws(strsplit(if (!is.null(args$models)) args$models else models_env, ",", fixed = TRUE)[[1]])
  models_sel <- models_sel[nzchar(models_sel)]

  # Read data and align numeric columns
  Xtr_all <- read_numeric_matrix(train_csv, n_rows = n_train)
  Xte_all <- read_numeric_matrix(test_csv,  n_rows = n_test)
  cn_tr <- colnames(Xtr_all); cn_te <- colnames(Xte_all)
  common <- intersect(cn_tr, cn_te)
  if (!length(common)) stop("No common numeric columns between train/test CSVs")
  X_tr <- Xtr_all[, common, drop = FALSE]
  X_te <- Xte_all[, common, drop = FALSE]
  stopifnot(is.matrix(X_tr), is.matrix(X_te), ncol(X_tr) == ncol(X_te))
  K <- ncol(X_tr)

  S <- list(X_tr = X_tr, X_te = X_te)
  cfg <- replicate(K, list(distr = "norm"), simplify = FALSE)

  # Composite nutzt interne Defaults/Val und ggf. Tuning

  # Train models with timing
  t_true_tr <- t_trtf_tr <- t_ttm_tr <- t_sep_tr <- t_cop_tr <- NA_real_
  t_true_te <- t_trtf_te <- t_ttm_te <- t_sep_te <- t_cop_te <- NA_real_

  mods <- list()

  if (want_model("true", models_sel)) {
    t_true_tr <- system.time(mods$true <- fit_TRUE(S, cfg))[['elapsed']]
    t_true_te <- system.time({ invisible(logL_TRUE(mods$true, X_te)) })[['elapsed']]
  }
  if (want_model("trtf", models_sel)) {
    mods$trtf <- tryCatch({
      t0 <- system.time({ model_trtf <- fit_TRTF(S, cfg, seed = seed) })
      t_trtf_tr <<- t0[['elapsed']]
      model_trtf
    }, error = function(e) {
      message("[WARN] TRTF unavailable: ", e$message)
      NULL
    })
    if (!is.null(mods$trtf)) t_trtf_te <- tryCatch(system.time(predict(mods$trtf, X_te, type = "logdensity_by_dim"))[["elapsed"]], error = function(e) NA_real_)
  }
  if (want_model("ttm", models_sel)) {
    fit_m <- fit_ttm(S, algo = "marginal", seed = seed); mods$ttm <- fit_m$S; t_ttm_tr <- fit_m$time_train
    t_ttm_te <- system.time(predict_ttm(mods$ttm, X_te, type = "logdensity_by_dim"))[["elapsed"]]
  }
  if (want_model("ttm_sep", models_sel)) {
    fit_s <- fit_ttm(S, algo = "separable", seed = seed); mods$ttm_sep <- fit_s$S; t_sep_tr <- fit_s$time_train
    t_sep_te <- system.time(predict_ttm(mods$ttm_sep, X_te, type = "logdensity_by_dim"))[["elapsed"]]
  }
  # Composite entfernt (nur TRTF, TTM marginal, TTM separable, Copula NP)
  if (want_model("copula_np", models_sel) && exists("fit_copula_np")) {
    t_cop_tr <- tryCatch(system.time(mods$copula_np <- fit_copula_np(S, seed = seed))[["elapsed"]], error = function(e) NA_real_)
    t_cop_te <- tryCatch(system.time(predict(mods$copula_np, X_te, type = "logdensity_by_dim"))[["elapsed"]], error = function(e) NA_real_)
  }

  # Build evaluation table (per-dimension NLL with SEs; sum-row has joint SE)
  N <- nrow(X_te)
  dims <- as.character(seq_len(K))
  fmt <- function(m, se) sprintf("%.2f ± %.2f", round(m, 2), round(2 * se, 2))

  tab <- data.frame(dim = dims, distribution = rep("data", K), stringsAsFactors = FALSE)

  # TRUE (marginal) baseline
  if (!is.null(mods$true)) {
    LD_true <- do.call(cbind, lapply(seq_len(K), function(k) .log_density_vec(X_te[, k], cfg[[k]]$distr, mods$true$theta[[k]])))
    ll_true <- -LD_true
    tab[["True (marginal)"]] <- fmt(colMeans(ll_true), apply(ll_true, 2, stderr))
    se_sum_true <- stats::sd(rowSums(ll_true)) / sqrt(N)
  } else {
    se_sum_true <- NA_real_
  }

  # TRTF
  if (!is.null(mods$trtf)) {
    ll_trtf <- tryCatch(-predict(mods$trtf, X_te, type = "logdensity_by_dim"), error = function(e) matrix(NA_real_, N, K))
    tab[["Random Forest"]] <- fmt(colMeans(ll_trtf), apply(ll_trtf, 2, stderr))
    se_sum_trtf <- stats::sd(rowSums(ll_trtf)) / sqrt(N)
  } else {
    se_sum_trtf <- NA_real_
  }

  # TTM variants
  if (!is.null(mods$ttm)) {
    ll_ttm <- -predict_ttm(mods$ttm, X_te, type = "logdensity_by_dim")
    tab[["Marginal Map"]] <- fmt(colMeans(ll_ttm), apply(ll_ttm, 2, stderr))
    se_sum_ttm <- stats::sd(rowSums(ll_ttm)) / sqrt(N)
  } else { se_sum_ttm <- NA_real_ }

  if (!is.null(mods$ttm_sep)) {
    ll_sep <- -predict_ttm(mods$ttm_sep, X_te, type = "logdensity_by_dim")
    tab[["Separable Map"]] <- fmt(colMeans(ll_sep), apply(ll_sep, 2, stderr))
    se_sum_sep <- stats::sd(rowSums(ll_sep)) / sqrt(N)
  } else { se_sum_sep <- NA_real_ }


  if (!is.null(mods$copula_np)) {
    ll_cop <- tryCatch(-predict(mods$copula_np, X_te, type = "logdensity_by_dim"), error = function(e) matrix(NA_real_, N, K))
    tab[["Copula NP"]] <- fmt(colMeans(ll_cop), apply(ll_cop, 2, stderr))
    se_sum_cop <- stats::sd(rowSums(ll_cop)) / sqrt(N)
  } else { se_sum_cop <- NA_real_ }

  # Sum row with joint SE for each model
  sum_row <- data.frame(dim = "k", distribution = "SUM", stringsAsFactors = FALSE)
  if (!is.null(mods$true))     sum_row[["True (marginal)"]] <- fmt(sum(colMeans(ll_true)), se_sum_true)
  if (!is.null(mods$trtf))     sum_row[["Random Forest"]]   <- fmt(sum(colMeans(ll_trtf)), se_sum_trtf)
  if (!is.null(mods$ttm))      sum_row[["Marginal Map"]]    <- fmt(sum(colMeans(ll_ttm)), se_sum_ttm)
  if (!is.null(mods$ttm_sep))  sum_row[["Separable Map"]]   <- fmt(sum(colMeans(ll_sep)),  se_sum_sep)
  if (!is.null(mods$copula_np))sum_row[["Copula NP"]]       <- fmt(sum(colMeans(ll_cop)),  se_sum_cop)
  tab <- rbind(tab, sum_row)

  # Output
  cat(sprintf("Train CSV: %s\nTest CSV: %s\n", train_csv, test_csv))
  cat(sprintf("n_train=%d, n_test=%d, K=%d\n", nrow(X_tr), nrow(X_te), K))
  print(tab)

  # Timing table mirroring main.R (unavailable entries kept NA)
  time_tab <- data.frame(
    model = c("True (marginal)", "Random Forest", "Copula NP", "Marginal Map", "Separable Map"),
    train_sec = c(t_true_tr, t_trtf_tr, t_cop_tr, t_ttm_tr, t_sep_tr),
    test_sec  = c(t_true_te, t_trtf_te, t_cop_te, t_ttm_te, t_sep_te),
    stringsAsFactors = FALSE
  )
  time_tab$total_sec <- with(time_tab, train_sec + test_sec)
  print(time_tab)
  timing_table <<- time_tab
  results_table <<- tab
  invisible(tab)
}

if (sys.nframe() == 0L) {
  invisible(main_power())
}

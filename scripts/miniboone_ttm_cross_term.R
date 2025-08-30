#!/usr/bin/env Rscript
# Grid over TTM Cross-term on MiniBooNE CSVs (small n)
# - User can set k (first columns) and n (first rows) right below
# - Q ∈ {6K, 8K, 10K}
# - λ_non,0 = 0.05 * n / Q,  λ_mon,0 = 2^-1 * λ_non,0
# - Reports table with hyperparameters, nats (SUM NLL on validation), and runtime; ranked by nats asc

# --- User-editable selection (set here) ---
k <- 4L   # keep first k columns; set to NA to keep all
n <- 50L  # use first n rows from each split
# -----------------------------------------

options(error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 1) })

suppressMessages({
  source("00_globals.R")
  source("models/ttm_separable.R")
  source("models/ttm_cross_term.R")
})

# Helper: read first N numeric rows from CSV
read_first_n_numeric <- function(path, n) {
  if (!file.exists(path)) stop("CSV not found: ", path)
  df <- read.csv(path, check.names = FALSE)
  keep <- vapply(df, is.numeric, logical(1))
  X <- as.matrix(df[seq_len(min(n, nrow(df))), keep, drop = FALSE])
  storage.mode(X) <- "double"
  if (!is.matrix(X) || any(!is.finite(X))) stop("Non-numeric or non-finite entries in ", path)
  X
}

stderr <- function(x) stats::sd(x)/sqrt(length(x))
sum_nll <- function(LD_by_dim) sum(colMeans(-LD_by_dim))

main <- function() {
  set.seed(as.integer(Sys.getenv("SEED", "42")))
  # Use n from header; fallback to 50 if invalid
  n_use <- suppressWarnings(as.integer(n))
  if (is.na(n_use) || !is.finite(n_use) || n_use <= 0L) n_use <- 50L
  # Use k from header; NA means keep all
  D_KEEP <- suppressWarnings(as.integer(k))
  train_csv <- Sys.getenv("TRAIN_CSV", "data/miniboone_train.csv")
  val_csv   <- Sys.getenv("VAL_CSV",   "data/miniboone_val.csv")
  test_csv  <- Sys.getenv("TEST_CSV",  "data/miniboone_test.csv")

  Xtr <- read_first_n_numeric(train_csv, n_use)
  Xval <- read_first_n_numeric(val_csv,   n_use)
  Xte <- read_first_n_numeric(test_csv,   n_use)
  stopifnot(ncol(Xtr) == ncol(Xval), ncol(Xtr) == ncol(Xte))
  if (!is.na(D_KEEP) && is.finite(D_KEEP) && D_KEEP >= 1L && D_KEEP < ncol(Xtr)) {
    cols <- seq_len(D_KEEP)
    Xtr  <- Xtr[, cols, drop = FALSE]
    Xval <- Xval[, cols, drop = FALSE]
    Xte  <- Xte[, cols, drop = FALSE]
  }
  # Train-only per-column log transform: y = log(x + c_k), c_k ensures positivity on train
  # Compute shift from train only and apply to all splits (no leakage of values)
  train_min <- apply(Xtr, 2, min)
  shift <- ifelse(train_min > 0, 0, -train_min + 1e-6)
  add_shift <- function(X, s) sweep(X, 2, s, "+")
  eps <- 1e-12
  Xtr <- log(pmax(add_shift(Xtr, shift), eps))
  Xval <- log(pmax(add_shift(Xval, shift), eps))
  Xte <- log(pmax(add_shift(Xte, shift), eps))
  n <- nrow(Xtr); K <- ncol(Xtr)
  Qs <- as.integer(c(6L*K, 8L*K, 10L*K))
  multipliers <- c(0.5, 1, 2)

  dir.create("results", showWarnings = FALSE)
  rows <- list()
  best <- list(nats = Inf, Q = NA_integer_, lambda_non = NA_real_, lambda_mon = NA_real_, idx = NA_integer_)

  for (Q in Qs) {
    # Regularisierung je Q: λ_non,0 = 0.05 * Q / n; λ_mon,0 = 0.5 * λ_non,0
    base_non <- 0.05 * Q / n
    base_mon <- 0.5 * base_non
    for (m in multipliers) {
      lambda_non0 <- base_non * m
      lambda_mon0 <- base_mon * m
      S <- list(X_tr = Xtr, X_te = Xval)
      t0 <- Sys.time()
      fit <- tryCatch(trainCrossTermMap(S, degree_g = 3, Q = Q,
                                        lambda_non = lambda_non0, lambda_mon = lambda_mon0,
                                        warmstart_from_separable = TRUE,
                                        sep_degree_g = 2, sep_lambda = 1e-3,
                                        seed = 42),
                      error = function(e) e)
      runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (inherits(fit, "error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          n = n, K = K, Q = Q, mult = m,
          lambda_non = lambda_non0, lambda_mon = lambda_mon0,
          nats = NA_real_, runtime_sec = runtime,
          stringsAsFactors = FALSE
        )
        next
      }
      LD_val <- tryCatch(-predict(fit$S, Xval, type = "logdensity_by_dim"), error = function(e) e)
      if (inherits(LD_val, "error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          n = n, K = K, Q = Q, mult = m,
          lambda_non = lambda_non0, lambda_mon = lambda_mon0,
          nats = NA_real_, runtime_sec = runtime,
          stringsAsFactors = FALSE
        )
        next
      }
      stopifnot(is.matrix(LD_val), nrow(LD_val) == nrow(Xval), ncol(LD_val) == K)
      nats <- sum(colMeans(LD_val))
      rows[[length(rows) + 1L]] <- data.frame(
        n = n, K = K, Q = Q, mult = m,
        lambda_non = lambda_non0, lambda_mon = lambda_mon0,
        nats = nats, runtime_sec = runtime,
        stringsAsFactors = FALSE
      )
      if (is.finite(nats) && nats < best$nats) best <- list(nats = nats, Q = Q, lambda_non = lambda_non0, lambda_mon = lambda_mon0, idx = length(rows))
    }
  }

  df <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()
  ord <- order(df$nats, df$runtime_sec, na.last = TRUE)
  df_ranked <- df[ord, c("n","K","Q","mult","lambda_non","lambda_mon","nats","runtime_sec") , drop = FALSE]
  rownames(df_ranked) <- NULL
  out_path <- sprintf("results/miniboone_ctm_grid_n%03d.csv", n)
  write.csv(df_ranked, out_path, row.names = FALSE)
  print(df_ranked)

  # Evaluate best on held-out test (optional)
  if (is.finite(best$nats)) {
    cat(sprintf("\n[Best] Q=%d, lambda_non=%.6g, lambda_mon=%.6g (val SUM NLL=%.4f)\n", best$Q, best$lambda_non, best$lambda_mon, best$nats))
    # Reuse best fit if index tracked; otherwise refit
    S_test <- list(X_tr = Xtr, X_te = Xte)
    fit_best <- tryCatch(trainCrossTermMap(S_test, degree_g = 3, Q = best$Q,
                                           lambda_non = best$lambda_non, lambda_mon = best$lambda_mon,
                                           warmstart_from_separable = TRUE,
                                           sep_degree_g = 2, sep_lambda = 1e-3,
                                           seed = 42), error = function(e) e)
    if (!inherits(fit_best, "error")) {
      LD_test <- -predict(fit_best$S, Xte, type = "logdensity_by_dim")
      nats_test <- sum(colMeans(LD_test))
      cat(sprintf("[Best→Test] SUM NLL (nats): %.4f\n", nats_test))
    } else {
      cat("[Best→Test] Fit failed: ", conditionMessage(fit_best), "\n", sep = "")
    }
  }
}

if (sys.nframe() == 0L) main()

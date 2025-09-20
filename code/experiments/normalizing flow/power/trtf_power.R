### Minimaler, schneller Konfigurationsblock (oben anpassbar) #################
# A) Anzahl Zeilen pro Split (Train/Test)
N_ROWS <- 10000L
# B) CSV-Pfade relativ zum Skript-Verzeichnis
SCRIPT_DIR <- (function(){
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))))
  # when sourced
  f <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
  if (is.character(f) && length(f) == 1L && !is.na(f)) return(dirname(f))
  getwd()
})()

locate_data_dir <- function(initial_dir, required) {
  candidates <- unique(c(
    initial_dir,
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "code", "experiments", "normalizing flow", "power"), winslash = "/", mustWork = FALSE)
  ))
  for (dir in candidates) {
    if (nzchar(dir) && all(file.exists(file.path(dir, required)))) return(dir)
  }
  initial_dir
}

SCRIPT_DIR <- locate_data_dir(SCRIPT_DIR, c("power_train_10000.csv", "power_test_10000.csv"))
TRAIN_CSV <- file.path(SCRIPT_DIR, "power_train_10000.csv")
TEST_CSV  <- file.path(SCRIPT_DIR, "power_test_10000.csv")
# C) Pfad für Evaluationstext (Konsole + Datei)
RESULTS_TXT <- file.path(SCRIPT_DIR, "results", "trtf_power_eval.txt")
# D) Hyperparameter
HP <- list(
  ntree     = 200L,
  minsplit  = 100L,
  minbucket = 20L,
  maxdepth  = 4L,
  seed      = 42L
)
###############################################################################

ensure_required_packages <- function(pkgs) {
  missing <- character(0)
  for (pkg in pkgs) {
    ok <- suppressPackageStartupMessages(require(pkg, quietly = TRUE, character.only = TRUE))
    if (!isTRUE(ok)) missing <- c(missing, pkg)
  }
  if (length(missing)) {
    warning(sprintf('Missing required packages: %s (continuing, but functionality may fail).',
                    paste(missing, collapse = ', ')))
    return(invisible(FALSE))
  }
}

ensure_required_packages(c('tram', 'trtf', 'partykit', 'mlt'))

mclapply_trtf <- function(X, FUN, cores = 1L, ...) {
  parallel::mclapply(X, FUN, ..., mc.cores = cores)
}

# Dynamische Parallelisierung (alle Kerne minus zwei)
detect_cores_safe <- function() {
  if (requireNamespace("parallel", quietly = TRUE)) {
    val <- tryCatch(parallel::detectCores(), error = function(e) 1L)
    if (!is.finite(val) || is.na(val) || val < 1L) 1L else as.integer(val)
  } else 1L
}
NC <- { tot <- detect_cores_safe(); max(1L, tot - 2L ) }
# Globale mc.cores-Option nicht setzen; Steuerung erfolgt ueber get_train_cores()
# options(mc.cores = NC)
NC
# Einheitliche Steuerung der Trainingskerne
get_train_cores <- function(default_max = 8L) {
  req <- as.integer(getOption("trtf.train_cores", default_max))
  dc  <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  if (!is.finite(dc) || is.na(dc) || dc < 1L) dc <- 1L
  req <- if (!is.finite(req) || is.na(req) || req < 1L) default_max else req
  min(req, as.integer(dc))
}

# Konfigurierbare Optionen (Defaults per getOption):
# - trtf.chunk_size:   Chunkgröße für Vorhersage (Default 5000L)
# - trtf.cores_fit_max: Max. Trainingskerne für traforest (Default 2L)
# - trtf.use_cli:      CLI-Progress nutzen (Default: Autodetect via requireNamespace)
# - trtf.rss_limit_mb: RSS-Grenze in MB für Memory-Guard (Default 12288L)
# - trtf.varimp:       VarImp berechnen (Default FALSE; opt-in)

# Optional: hübschere Fortschrittsbalken via 'cli' falls vorhanden
{ .cli_avail <- requireNamespace("cli", quietly = TRUE)
  .use_cli <- isTRUE(getOption("trtf.use_cli", .cli_avail)) && .cli_avail }

# Begrenze BLAS/OMP Threads auf 1 (falls nicht gesetzt)
if (!nzchar(Sys.getenv("OMP_NUM_THREADS", "")))          Sys.setenv(OMP_NUM_THREADS = "1")
if (!nzchar(Sys.getenv("OPENBLAS_NUM_THREADS", "")))      Sys.setenv(OPENBLAS_NUM_THREADS = "1")
if (!nzchar(Sys.getenv("MKL_NUM_THREADS", "")))          Sys.setenv(MKL_NUM_THREADS = "1")
if (!nzchar(Sys.getenv("VECLIB_MAXIMUM_THREADS", "")))   Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")

# Robuste Fork-Erkennung (ohne direkte Abhängigkeit auf supportsMulticore Symbol)
has_fork <- function() {
  if (!requireNamespace("parallel", quietly = TRUE)) return(FALSE)
  can_fork <- tryCatch({
    f <- get("supportsMulticore", asNamespace("parallel"))
    isTRUE(f())
  }, error = function(e) .Platform$OS.type == "unix")
  can_fork && Sys.getenv("RSTUDIO") != "1"
}

# Inline TRTF Modellcode (ohne weitere Abhängigkeiten aus diesem Repo)
trtf_params <- list(minsplit = HP$minsplit, minbucket = HP$minbucket, maxdepth = HP$maxdepth, seed = HP$seed)

mytrtf <- function(data, ntree, minsplit, minbucket, maxdepth, seed, cores = NC,
                   verbose = TRUE, show_eta = TRUE, trace_trees = FALSE) {
  if (!is.matrix(data)) {
    warning("Input 'data' is not a matrix; attempting to coerce.")
    data <- as.matrix(data)
  }
  set.seed(seed)
  K <- ncol(data)
  df <- as.data.frame(data)
  names(df) <- paste0("X", seq_len(K))

  ymod <- lapply(names(df), function(y) {
    tram::BoxCox(as.formula(paste(y, "~ 1")), data = df)
  })

  forests <- vector("list", K - 1L)

  ctrl <- partykit::ctree_control(
    minsplit  = minsplit,
    minbucket = minbucket,
    maxdepth  = maxdepth,
    saveinfo  = FALSE
  )

  # Backend detection and training cores (fork vs PSOCK)
  supports_mc <- has_fork()
  cores_fit <- get_train_cores(8L)

  # --- Fortschritt (coarse) ---
  if (verbose) {
    msg <- sprintf("Fitting %d conditional forests (ntree=%d, depth<=%d) ...",
                   K - 1L, as.integer(ntree), as.integer(maxdepth))
    message(msg)
    if (.use_cli) {
      pb_id <- cli::cli_progress_bar(
        "Fitting conditional forests",
        total = K - 1L
      )
    } else {
      pb <- utils::txtProgressBar(min = 0, max = K - 1L, style = 3)
    }
  }
  t_start <- proc.time()[3]

  # cores_fit already determined above via get_train_cores()/detectCores fallback

  # Memory-Guard Setup (optional ps)
  .ps_avail <- requireNamespace("ps", quietly = TRUE)
  rss_peak_mb <- NA_real_
  rss_limit_mb <- {
    rl <- suppressWarnings(as.numeric(getOption("trtf.rss_limit_mb", 12288L)))
    if (!is.finite(rl) || is.na(rl) || rl <= 0) 12288 else rl
  }

  if (supports_mc) {
    # Forking available: parallelize inside each forest via mclapply
    for (k in 2:K) {
      rhs <- paste(names(df)[1:(k - 1)], collapse = "+")
      fm <- as.formula(paste(names(df)[k], "~", rhs))
      current_mtry <- max(1, floor((k - 1) / 2))

      t0 <- proc.time()[3]
      forests[[k - 1L]] <- trtf::traforest(
        ymod[[k]], formula = fm, data = df,
        trace = isTRUE(trace_trees), ntree = ntree,
        mltargs = list(), mtry = current_mtry,
        control = ctrl,
        applyfun = mclapply_trtf, cores = cores_fit
      )

      # Memory-Guard nach jedem Forest (optional, falls 'ps' vorhanden)
      if (.ps_avail) {
        rss_now <- try({
          as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]) / (1024^2)
        }, silent = TRUE)
        if (is.numeric(rss_now) && length(rss_now) == 1L && is.finite(rss_now)) {
          rss_peak_mb <- if (is.na(rss_peak_mb)) rss_now else max(rss_peak_mb, rss_now)
          if (rss_now > rss_limit_mb) {
            warning(sprintf("RSS limit exceeded: current=%.0f MB > limit=%.0f MB. Continuing anyway.",
                            rss_now, rss_limit_mb))
          }
        }
      }
      if (verbose) {
        done <- k - 1L
        if (.use_cli) {
          if (show_eta && done >= 1L) {
            elapsed <- proc.time()[3] - t_start
            rate    <- elapsed / done
            remain  <- rate * ((K - 1L) - done)
            cli::cli_progress_update(
              id = pb_id, set = done,
              status = sprintf("last %.1fs | ETA ~ %.1fs", proc.time()[3] - t0, remain)
            )
          } else {
            cli::cli_progress_update(id = pb_id, set = done)
          }
        } else {
          utils::setTxtProgressBar(pb, done)
          if (show_eta && done >= 1L) {
            elapsed <- proc.time()[3] - t_start
            rate    <- elapsed / done
            remain  <- rate * ((K - 1L) - done)
            cat(sprintf("  [Forest %d/%d: %.1fs | ETA ~ %.1fs]\n",
                        done, K - 1L, proc.time()[3] - t0, remain))
            flush.console()
          }
        }
      }
    }
  } else {
    # No forking (e.g., RStudio): parallelize across forests via PSOCK cluster
    workers <- {
      dc <- tryCatch(parallel::detectCores(), error = function(e) 1L)
      if (!is.finite(dc) || is.na(dc) || dc < 1L) dc <- 1L
      min(cores_fit, K - 1L, as.integer(dc))
    }
    cl <- parallel::makeCluster(workers, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { library(tram); library(trtf); library(partykit); library(mlt) })
    parallel::clusterExport(cl, c("df", "ymod", "ntree", "ctrl"), envir = environment())

    if (requireNamespace("pbapply", quietly = TRUE)) {
      forests_list <- pbapply::pblapply(
        2:K,
        function(k, df, ymod, ntree, ctrl) {
          rhs <- paste(names(df)[1:(k - 1)], collapse = "+")
          fm  <- as.formula(paste(names(df)[k], "~", rhs))
          current_mtry <- max(1, floor((k - 1) / 2))
          trtf::traforest(
            ymod[[k]], formula = fm, data = df,
            ntree = ntree, mtry = current_mtry, control = ctrl,
            trace = FALSE, mltargs = list(),
            applyfun = lapply, cores = 1L
          )
        },
        cl = cl
      )
    } else {
      forests_list <- parallel::parLapply(
        cl, 2:K,
        function(k, df, ymod, ntree, ctrl) {
          rhs <- paste(names(df)[1:(k - 1)], collapse = "+")
          fm  <- as.formula(paste(names(df)[k], "~", rhs))
          current_mtry <- max(1, floor((k - 1) / 2))
          trtf::traforest(
            ymod[[k]], formula = fm, data = df,
            ntree = ntree, mtry = current_mtry, control = ctrl,
            trace = FALSE, mltargs = list(),
            applyfun = lapply, cores = 1L
          )
        },
        df, ymod, ntree, ctrl
      )
    }

    forests <- vector("list", K - 1L)
    for (i in seq_along(forests_list)) forests[[i]] <- forests_list[[i]]

    # Mark progress as done
    if (verbose) {
      if (.use_cli) {
        cli::cli_progress_update(id = pb_id, set = K - 1L)
      } else {
        utils::setTxtProgressBar(pb, K - 1L)
      }
    }
  }
  if (verbose) {
    if (.use_cli) cli::cli_progress_done(id = pb_id) else close(pb)
  }

  # Backend summary
  message(sprintf("Backend: %s | train_cores=%d",
                  if (supports_mc) "mclapply" else "PSOCK",
                  as.integer(cores_fit)))

  # VarImp nur optional (opt-in) und nach Moeglichkeit subsampled
  vimp <- NULL
  if (isTRUE(getOption("trtf.varimp", FALSE))) {
    vimp <- lapply(forests, function(fr) {
      if (requireNamespace("partykit", quietly = TRUE)) {
        # Falls moeglich, subsample Daten fuer VarImp (ohne harte Abhaengigkeit)
        vi <- try({
          vfun <- get("varimp", asNamespace("partykit"))
          # Bestmoeglich ohne Annahmen ueber Methodensignatur
          do.call(vfun, list(fr))
        }, silent = TRUE)
        vi
      } else NULL
    })
  }

  # Ressourcen-Zusammenfassung
  if (verbose) {
    cs <- suppressWarnings(as.integer(getOption("trtf.chunk_size", 5000L)))
    if (!is.finite(cs) || is.na(cs) || cs < 1L) cs <- 5000L
    msg_sum <- sprintf("Resources: cores_fit=%d | chunk_size=%d%s",
                       as.integer(cores_fit), as.integer(cs),
                       if (!is.na(rss_peak_mb)) sprintf(" | rss_peak=%.0f MB", rss_peak_mb) else "")
    message(msg_sum)
  }

  res <- list(ymod = ymod, forests = forests, seed = seed, varimp = vimp)
  class(res) <- "mytrtf"
  res
}

predict.mytrtf <- function(object, newdata,
                           type = c("logdensity", "logdensity_by_dim"),
                           cores = NC, trace = FALSE) {
  type <- match.arg(type)
  if (!inherits(object, "mytrtf")) {
    warning("Predict called with non-'mytrtf' object; proceeding but results may be invalid.")
  }
  if (!is.matrix(newdata)) {
    warning("New data is not a matrix; attempting to coerce.")
    newdata <- as.matrix(newdata)
  }
  K <- length(object$ymod)
  # Standardize newdata if model carries a standardizer; set Jacobian per-dimension
  if (!is.null(object$std) && is.list(object$std)) {
    Xs <- apply_standardizer(newdata, object$std)
    jac_sigma <- -log(object$std$sd)
  } else {
    Xs <- newdata
    jac_sigma <- rep(0, K)
  }
  dfX <- as.data.frame(Xs)
  names(dfX) <- paste0("X", seq_len(K))
  N <- nrow(dfX)

  # Helper: chunked evaluation to avoid NxN allocations
  chunked_ld_j <- function(fr, df_all, resp, chunk_size = 5000L) {
    # chunk size via option override
    cs_opt <- getOption("trtf.chunk_size", chunk_size)
    cs <- as.integer(cs_opt)
    if (!is.finite(cs) || is.na(cs) || cs < 1L) cs <- as.integer(chunk_size)
    total <- N
    out <- numeric(total)
    n_chunks <- as.integer(ceiling(total / cs))
    show_progress <- isTRUE(trace)

    # progress setup per dimension
    if (show_progress) {
      if (.use_cli) {
        pbid <- cli::cli_progress_bar(
          "Predicting logdensity (chunks)", total = n_chunks
        )
      } else {
        pb <- utils::txtProgressBar(min = 0, max = n_chunks, style = 3)
      }
    }

    i <- 1L; cidx <- 0L
    while (i <= total) {
      j <- min(i + cs - 1L, total)
      idx <- i:j
      # vapply across rows in the chunk; predict with per-row data and scalar q
      t0 <- proc.time()[3]
      out[idx] <- vapply(idx, function(ii) {
        qii <- df_all[[resp]][ii]
        pii <- predict(fr,
                       newdata = df_all[ii, , drop = FALSE],
                       type = "logdensity", q = qii,
                       cores = 1L, trace = FALSE)
        as.numeric(pii)
      }, numeric(1))

      cidx <- cidx + 1L
      if (show_progress) {
        if (.use_cli) {
          # show last chunk duration and rough ETA
          elapsed_chunk <- proc.time()[3] - t0
          remain <- (n_chunks - cidx) * elapsed_chunk
          cli::cli_progress_update(id = pbid, set = cidx,
                                   status = sprintf("last %.1fs | ETA ~ %.1fs",
                                                    elapsed_chunk, remain))
        } else {
          utils::setTxtProgressBar(pb, cidx)
        }
      }
      i <- j + 1L
    }
    if (show_progress) {
      if (.use_cli) cli::cli_progress_done(id = pbid) else close(pb)
    }
    out
  }

  ld1 <- predict(object$ymod[[1]], newdata = dfX, type = "logdensity")
  ld1 <- as.numeric(ld1) + jac_sigma[1L]
  if (!is.numeric(ld1) || length(ld1) != N) {
    warning("Primary marginal logdensity has unexpected shape; coercing to numeric vector with length N.")
    ld1 <- rep(as.numeric(ld1), length.out = N)
  }

  # Parallel lapply over forests (fork if available, else serial)
  plapply_forests <- function(X, FUN, cores) {
    if (has_fork()) parallel::mclapply(X, FUN, mc.cores = min(cores, length(X)))
    else lapply(X, FUN)
  }

  # mclapply across forests; within each forest remain serial + chunked (no NxN)
  ld_rest <- plapply_forests(seq_along(object$forests), function(j) {
    fr <- object$forests[[j]]
    resp <- paste0("X", j + 1L)
    out <- chunked_ld_j(fr, dfX, resp, chunk_size = getOption("trtf.chunk_size", 5000L))
    out + jac_sigma[j + 1L]
  }, cores = if (exists("get_train_cores", inherits = TRUE)) get_train_cores(8L) else {
    dc <- tryCatch(parallel::detectCores(), error = function(e) 1L)
    if (!is.finite(dc) || is.na(dc) || dc < 1L) dc <- 1L
    min(8L, as.integer(dc))
  })

  ld_rest <- do.call(cbind, ld_rest)
  ll <- cbind(ld1, ld_rest)

  if (!is.matrix(ll) || nrow(ll) != N || ncol(ll) != K) {
    warning("Per-dimension logdensity output has unexpected shape; attempting to coerce to matrix.")
    ll <- matrix(as.numeric(ll), nrow = N, ncol = K)
  }
  if (!all(is.finite(ll))) {
    warning("Non-finite values in TRTF logdensity_by_dim; replacing with NA.")
    ll[!is.finite(ll)] <- NA_real_
  }
  if (type == "logdensity_by_dim") return(ll)
  joint <- rowSums(ll)
  if (!all(is.finite(joint))) {
    warning("Non-finite values in TRTF joint logdensity; replacing with NA.")
    joint[!is.finite(joint)] <- NA_real_
  }
  check_joint_logdensity(joint, ll, warn_abs = 1e-10, warn_rel = 1e-8)
  joint
}

logL_TRTF <- function(model, X, cores = NC) {
  nll_value <- -mean(predict(model, X, type = "logdensity",
                             cores = cores, trace = FALSE))
  if (!is.finite(nll_value)) {
    warning("log-likelihood not finite; returning NA.")
    return(NA_real_)
  }
  nll_value
}

logL_TRTF_dim <- function(model, X, cores = NC) {
  ll <- predict(model, X, type = "logdensity_by_dim",
                cores = cores, trace = FALSE)
  res <- -colMeans(ll)
  if (!all(is.finite(res))) {
    warning("Per-dimension log-likelihood not finite; returning NA.")
    res[!is.finite(res)] <- NA_real_
  }
  res
}

fit_TRTF <- function(S, config, seed = NULL, cores = NC) {
  if (!is.list(S)) {
    warning("Input S is not a list; attempting to wrap into list structure.")
    S <- as.list(S)
  }
  X_tr <- S$X_tr
  X_te <- S$X_te
  if (!is.matrix(X_tr)) {
    warning("Training data is not a matrix; attempting to coerce.")
    X_tr <- as.matrix(X_tr)
  }
  if (!is.matrix(X_te)) {
    warning("Test data is not a matrix; attempting to coerce.")
    X_te <- as.matrix(X_te)
  }

  if (!is.null(seed)) set.seed(seed)

  # train-only standardization
  std   <- make_standardizer(X_tr)
  X_trs <- apply_standardizer(X_tr, std)
  X_tes <- apply_standardizer(X_te, std)

  # ntree aus HP; train on standardized data
  mod <- mytrtf(data = X_trs,
                ntree = as.integer(HP$ntree),
                minsplit = trtf_params$minsplit,
                minbucket = trtf_params$minbucket,
                maxdepth = trtf_params$maxdepth,
                seed = seed,
                cores = cores)
  
  mod$config  <- config
  mod$std     <- std

  # --- Metrics: store LL and NLL consistently ---
  LL_train      <- mean(predict(mod, X_tr, type = "logdensity", cores = cores, trace = FALSE))
  LL_te_by_dim  <- predict(mod, X_te, type = "logdensity_by_dim", cores = cores, trace = FALSE)
  LL_test       <- sum(colMeans(LL_te_by_dim))
  SE_LL_test    <- stats::sd(rowSums(LL_te_by_dim)) / sqrt(nrow(X_te))

  mod$LL_train    <- LL_train
  mod$NLL_train   <- -LL_train
  mod$LL_test     <- LL_test
  mod$NLL_test    <- -LL_test
  mod$SE_LL_test  <- SE_LL_test
  mod$SE_NLL_test <- SE_LL_test
  mod$LL_te_dim   <- colMeans(LL_te_by_dim)
  
  mod
}
# TRTF-only Runner (ohne TTM/Copula)

# Kleine Helfer
stderr <- function(x) stats::sd(x) / sqrt(length(x))

read_numeric_matrix <- function(path, n_rows = NA_integer_) {
  if (!file.exists(path)) {
    warning(sprintf("CSV not found: %s; returning empty matrix.", path))
    return(matrix(numeric(0), nrow = 0))
  }
  df <- utils::read.csv(path, header = TRUE, check.names = TRUE)
  # Keep only numeric columns present
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  if (!length(num_cols)) {
    warning(sprintf("No numeric columns in %s; returning empty matrix.", basename(path)))
    return(matrix(numeric(0), nrow = 0))
  }
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

# ---- train-only standardization ------------------------------------------------
make_standardizer <- function(Xtr) {
  mu <- colMeans(Xtr)
  sd <- apply(Xtr, 2, stats::sd)
  sd[!is.finite(sd) | is.na(sd) | sd <= 0] <- 1.0
  list(mu = as.numeric(mu), sd = as.numeric(sd))
}

apply_standardizer <- function(X, S) {
  Xs <- sweep(X, 2, S$mu, "-")
  Xs <- sweep(Xs, 2, S$sd, "/")
  storage.mode(Xs) <- "double"
  Xs
}

check_joint_logdensity <- function(joint, per_dim,
                                   warn_abs = 1e-8,
                                   warn_rel = 1e-6,
                                   stop_abs = 1e-4,
                                   stop_rel = 1e-4) {
  result <- list(max_abs = NA_real_, max_rel = NA_real_)
  if (!is.numeric(joint) || (!is.matrix(per_dim) && !is.data.frame(per_dim))) {
    warning("Joint/per-dimension logdensity check skipped: unexpected input types.")
    return(result)
  }
  per_dim <- as.matrix(per_dim)
  if (length(joint) != nrow(per_dim)) {
    warning("Joint/per-dimension logdensity check skipped: length mismatch.")
    return(result)
  }
  residual <- joint - rowSums(per_dim)
  if (!all(is.finite(residual))) {
    warning("Joint/per-dimension residuals are not finite; skipping strict consistency check.")
    return(result)
  }
  max_abs <- max(abs(residual))
  max_rel <- max(abs(residual) / pmax(1, abs(joint)))
  result$max_abs <- max_abs
  result$max_rel <- max_rel
  if (max_abs > warn_abs || max_rel > warn_rel) {
    msg <- sprintf(
      "Joint/per-dimension logdensity mismatch: max_abs=%.3e, max_rel=%.3e",
      max_abs, max_rel
    )
    if (max_abs > stop_abs && max_rel > stop_rel) {
      warning(paste0(msg, " (exceeds nominal tolerance, continuing)."))
    } else {
      message(msg)
    }
  }
  invisible(result)
}

trtf_power <- function() {
  seed <- as.integer(HP$seed)
  set.seed(seed)

  # CSV-Pfade aus dem Kopfbereich
  train_csv <- TRAIN_CSV
  test_csv  <- TEST_CSV

  # Anzahl Zeilen (Train/Test)
  n_train <- as.integer(N_ROWS)
  n_test  <- as.integer(N_ROWS)

  # Read data and align numeric columns
  Xtr_all <- read_numeric_matrix(train_csv, n_rows = n_train)
  Xte_all <- read_numeric_matrix(test_csv,  n_rows = n_test)
  cn_tr <- colnames(Xtr_all); cn_te <- colnames(Xte_all)
  common <- intersect(cn_tr, cn_te)
  if (!length(common)) {
    warning("No common numeric columns between train/test CSVs; using union of available columns with recycling.")
    max_cols <- max(length(cn_tr), length(cn_te))
    cn_tr <- rep(cn_tr, length.out = max_cols)
    cn_te <- rep(cn_te, length.out = max_cols)
    common <- unique(c(cn_tr, cn_te))
  }
  X_tr <- Xtr_all[, common, drop = FALSE]
  X_te <- Xte_all[, common, drop = FALSE]
  if (!is.matrix(X_tr)) {
    warning("X_tr is not a matrix after coercion; results may be invalid.")
  }
  if (!is.matrix(X_te)) {
    warning("X_te is not a matrix after coercion; results may be invalid.")
  }
  if (ncol(X_tr) != ncol(X_te)) {
    warning("Train/Test column counts differ; truncating to common columns.")
    common_cols <- seq_len(min(ncol(X_tr), ncol(X_te)))
    X_tr <- X_tr[, common_cols, drop = FALSE]
    X_te <- X_te[, common_cols, drop = FALSE]
  }
  K <- ncol(X_tr)

  S <- list(X_tr = X_tr, X_te = X_te)
  cfg <- replicate(K, list(distr = "norm"), simplify = FALSE)

  # Train models with timing
  t_trtf_tr <- t_trtf_te <- NA_real_

  cores_use <- NC
  mod_trtf <- tryCatch({
    t0 <- system.time({ m <- fit_TRTF(S, cfg, seed = seed, cores = cores_use) })
    t_trtf_tr <- t0[['elapsed']]
    m
  }, error = function(e) {
    warning(sprintf("[WARNING] TRTF unavailable: %s", e$message))
    return(NULL)
  })
  if (is.null(mod_trtf)) {
    warning("TRTF model could not be trained; downstream metrics will be set to NA.")
  }
  # --- Single test-predict (timed) without sign flip ---
  if (!is.null(mod_trtf)) {
    t_pred <- tryCatch(system.time({
      LL_by_dim <- predict(mod_trtf, X_te, type = "logdensity_by_dim",
                           cores = cores_use,
                           trace = isTRUE(getOption("trtf.trace_predict", FALSE)))
    }), error = function(e) {
      warning(sprintf("Prediction failed: %s", e$message))
      structure(NA_real_, names = c("user.self","sys.self","elapsed","user.child","sys.child"))
    })
    t_trtf_te <- unname(t_pred[["elapsed"]])
  } else {
    LL_by_dim <- matrix(NA_real_, nrow(X_te), ncol(X_te))
    t_trtf_te <- NA_real_
  }
  if (!exists("LL_by_dim", inherits = FALSE)) LL_by_dim <- matrix(NA_real_, nrow(X_te), ncol(X_te))

  # --- Invariants: joint equals row sums ---
  L_joint <- if (!is.null(mod_trtf)) {
    tryCatch(predict(mod_trtf, X_te, type = "logdensity", cores = cores_use, trace = FALSE),
             error = function(e) {
               warning(sprintf("Joint prediction failed: %s", e$message))
               rep(NA_real_, nrow(X_te))
             })
  } else {
    rep(NA_real_, nrow(X_te))
  }
  if (!is.matrix(LL_by_dim)) {
    warning("LL_by_dim not a matrix; coercing for downstream summaries.")
    LL_by_dim <- as.matrix(LL_by_dim)
  }
  if (length(L_joint) != nrow(LL_by_dim)) {
    warning("Joint predictions length does not match LL_by_dim rows; recycling values.")
    L_joint <- rep(L_joint, length.out = nrow(LL_by_dim))
  }
  joint_check <- check_joint_logdensity(L_joint, LL_by_dim)
  message(sprintf("Joint consistency check: max_abs=%.3e | max_rel=%.3e",
                  joint_check$max_abs, joint_check$max_rel))

  # --- Evaluation table: LL and NLL with 2*SE ---
  N <- nrow(X_te); K <- ncol(X_te)
  stderr <- function(x) stats::sd(x) / sqrt(length(x))
  fmt <- function(m, se) sprintf("%.2f +/- %.2f", round(m, 2), round(2 * se, 2))

  dims <- as.character(seq_len(K))
  tab <- data.frame(dim = dims, distribution = rep("data", K), stringsAsFactors = FALSE)

  LL_mean  <- colMeans(LL_by_dim)
  LL_se    <- apply(LL_by_dim, 2, stderr)
  tab[["LL (RF)"]]  <- fmt(LL_mean, LL_se)
  tab[["NLL (RF)"]] <- fmt(-LL_mean, LL_se)

  joint_LL <- mean(rowSums(LL_by_dim))
  joint_SE <- stats::sd(rowSums(LL_by_dim)) / sqrt(N)
  sum_row <- data.frame(dim = "k", distribution = "SUM", stringsAsFactors = FALSE)
  sum_row[["LL (RF)"]]  <- fmt(joint_LL, joint_SE)
  sum_row[["NLL (RF)"]] <- fmt(-joint_LL, joint_SE)
  tab <- rbind(tab, sum_row)

  # --- Print results and timing ---
  sink_path <- RESULTS_TXT
  sink_active <- FALSE
  if (!is.null(sink_path) && nzchar(sink_path)) {
    dir.create(dirname(sink_path), recursive = TRUE, showWarnings = FALSE)
    sink_attempt <- try(sink(sink_path, split = TRUE), silent = TRUE)
    if (inherits(sink_attempt, "try-error")) {
      cond <- attr(sink_attempt, "condition")
      warn_msg <- if (!is.null(cond)) conditionMessage(cond) else as.character(sink_attempt)
      warning(sprintf("Could not mirror output to %s: %s",
                      sink_path, warn_msg))
    } else {
      sink_active <- TRUE
      on.exit({ sink(NULL) }, add = TRUE)
    }
  }
  op <- options(width = max(200L, getOption("width")))
  on.exit(options(op), add = TRUE)
  print(tab, row.names = FALSE, right = TRUE)
  cat(sprintf(
    "hp: ntree=%d, minsplit=%d, minbucket=%d, maxdepth=%d, seed=%d, train_cores=%d\n",
    as.integer(HP$ntree), as.integer(HP$minsplit), as.integer(HP$minbucket),
    as.integer(HP$maxdepth), as.integer(HP$seed), get_train_cores(8L)
  ))

  # Timing table (TRTF only)
  time_tab <- data.frame(
    model = c("Random Forest"),
    train_sec = c(t_trtf_tr),
    test_sec  = c(t_trtf_te),
    stringsAsFactors = FALSE
  )
  time_tab$total_sec <- with(time_tab, train_sec + test_sec)
  print(time_tab)
  if (sink_active) {
    cat(sprintf("Output mirrored to %s\n",
                normalizePath(sink_path, winslash = "/", mustWork = FALSE)))
  }
  timing_table <<- time_tab
  results_table <<- tab
  invisible(tab)
}

if (sys.nframe() == 0L) {
  invisible(trtf_power())
}

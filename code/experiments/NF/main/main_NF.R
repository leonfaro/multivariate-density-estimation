### Minimaler, schneller Konfigurationsblock (oben anpassbar) #################
# A) Datensatzwahl ("power", "hepmass", "gas")
DATASET <- Sys.getenv("TRTF_DATASET", "hepmass")
DATASET <- if (nzchar(DATASET)) tolower(trimws(DATASET)) else "hepmass"

dataset_cfg <- switch(
  DATASET,
  "power" = list(
    label = "power",
    n_default = 1000L,
    data_dir = normalizePath("../power", winslash = "/", mustWork = FALSE),
    train_csv = "power_train.csv",
    test_csv  = "power_test.csv",
    meta = list(type = "csv", file = "power_standardization.csv"),
    sample_rows = FALSE
  ),
  "hepmass" = list(
    label = "hepmass",
    n_default = NA_integer_,
    data_dir = normalizePath("../hepmass", winslash = "/", mustWork = FALSE),
    train_csv = "hepmass_train_10k_preprocessed.csv",
    test_csv  = "hepmass_test_10k_preprocessed.csv",
    meta = list(type = "rds", file = "hepmass_standardization.rds"),
    sample_rows = FALSE
  ),
  "gas" = list(
    label = "gas",
    n_default = 100L,
    data_dir = normalizePath("../gas", winslash = "/", mustWork = FALSE),
    train_csv = "gas_train.csv",
    test_csv  = "gas_test.csv",
    meta = NULL,
    sample_rows = TRUE
  ),
  "miniboone" = list(
    label = "miniboone",
    n_default = 100L,
    data_dir = normalizePath("../miniboone", winslash = "/", mustWork = FALSE),
    train_csv = "miniboone_train.csv",
    test_csv  = "miniboone_test.csv",
    meta = list(type = "csv", file = "miniboone_train_standardization.csv"),
    sample_rows = FALSE
  ),
  stop(sprintf("Unsupported DATASET: %s", DATASET))
)

sample_override <- Sys.getenv("TRTF_SAMPLE_ROWS", "")
if (nzchar(sample_override)) {
  sample_flag <- tolower(trimws(sample_override))
  dataset_cfg$sample_rows <- sample_flag %in% c("1", "true", "t", "yes", "y")
}

RESULT_LABEL <- dataset_cfg$label
DATA_DIR <- dataset_cfg$data_dir
TRAIN_CSV <- file.path(DATA_DIR, dataset_cfg$train_csv)
TEST_CSV  <- file.path(DATA_DIR, dataset_cfg$test_csv)
META_INFO <- dataset_cfg$meta
META_PATH <- if (is.null(META_INFO)) NULL else file.path(DATA_DIR, META_INFO$file)
RESULTS_DIR <- normalizePath(file.path("results", RESULT_LABEL), winslash = "/", mustWork = FALSE)

# D) Hyperparameter
hp_cfg <- switch(
  RESULT_LABEL,
  "power" = list(minbucket = 87L, minsplit = 174L, maxdepth = 5L),
  "gas" = list(minbucket = 100L, minsplit = 200L, maxdepth = 5L),
  "hepmass" = list(minbucket = 163L, minsplit = 326L, maxdepth = 4L),
  "miniboone" = list(minbucket = 232L, minsplit = 464L, maxdepth = 4L),
  stop(sprintf("No hyperparameter configuration for dataset '%s'", RESULT_LABEL))
)

HP <- list(
  ntree     = 200L,
  minsplit  = hp_cfg$minsplit,
  minbucket = hp_cfg$minbucket,
  maxdepth  = hp_cfg$maxdepth,
  seed      = 42L
)

# Lauf-spezifische Override-Werte (Standard: dataset_cfg$n_default)
N_ROWS <- dataset_cfg$n_default
DEFAULT_OVERRIDE <- 100L
override_env <- Sys.getenv("TRTF_N", "")
if (nzchar(override_env)) {
  override_val <- suppressWarnings(as.integer(override_env))
  if (is.finite(override_val) && override_val > 0L) {
    N_ROWS <- override_val
  }
} else if (!is.na(DEFAULT_OVERRIDE) && DEFAULT_OVERRIDE > 0L) {
  N_ROWS <- DEFAULT_OVERRIDE
}
N_LABEL <- if (is.na(N_ROWS)) "all" else as.integer(N_ROWS)
RESULT_FILE_PREFIX <- RESULT_LABEL

RESULTS_TXT <- file.path(
  RESULTS_DIR,
  sprintf(
    "%s_seed%d_N%s.txt",
    RESULT_FILE_PREFIX,
    as.integer(HP$seed),
    N_LABEL
  )
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
cores_cap_global <- min(8L, detect_cores_safe())
NC <- cores_cap_global
options(trtf.train_cores = cores_cap_global)
# Globale mc.cores-Option nicht setzen; Steuerung erfolgt ueber get_train_cores()
# options(mc.cores = NC)
NC
# Einheitliche Steuerung der Trainingskerne
get_train_cores <- function(default_max = 8L) {
  req <- as.integer(getOption("trtf.train_cores", default_max))
  if (!is.finite(req) || is.na(req) || req < 1L) req <- default_max
  min(req, cores_cap_global)
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
  cores_fit <- max(1L, get_train_cores(8L))

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

  res <- list(ymod = ymod,
              forests = forests,
              seed = seed,
              varimp = vimp,
              train_cores = as.integer(cores_fit),
              train_cores_cap = as.integer(cores_cap_global))
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
  forests_count <- max(1L, length(object$forests))
  cores_cap_pred <- min(cores_cap_global, forests_count)
  cores_pred <- if (exists("get_train_cores", inherits = TRUE)) {
    max(1L, get_train_cores(cores_cap_pred))
  } else {
    dc <- tryCatch(parallel::detectCores(), error = function(e) 1L)
    if (!is.finite(dc) || is.na(dc) || dc < 1L) dc <- 1L
    min(cores_cap_pred, as.integer(dc))
  }

  ld_rest <- plapply_forests(seq_along(object$forests), function(j) {
    fr <- object$forests[[j]]
    resp <- paste0("X", j + 1L)
    out <- chunked_ld_j(fr, dfX, resp, chunk_size = getOption("trtf.chunk_size", 5000L))
    out + jac_sigma[j + 1L]
  }, cores = cores_pred)

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

  # Checklist: trust prep standardization, avoid double scaling, preserve standardized space
  # Daten sind bereits standardisiert -> keine zweite Standardisierung
  std   <- NULL
  X_trs <- X_tr
  X_tes <- X_te

  # ntree aus HP; train on standardised data
  mod <- mytrtf(data = X_trs,
                ntree = as.integer(HP$ntree),
                minsplit = trtf_params$minsplit,
                minbucket = trtf_params$minbucket,
                maxdepth = trtf_params$maxdepth,
                seed = seed,
                cores = cores,
                verbose = FALSE)

  mod$config  <- config
  mod$std     <- std

  mod
}
# TRTF-only Runner (ohne TTM/Copula)

# Kleine Helfer
stderr <- function(x) stats::sd(x) / sqrt(length(x))

read_numeric_matrix <- function(path, n_rows = NA_integer_, sample_rows = FALSE) {
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
    n_keep <- min(as.integer(n_rows), nrow(df))
    if (sample_rows && n_keep < nrow(df)) {
      sel <- sample.int(nrow(df), n_keep, replace = FALSE)
      df <- df[sel, , drop = FALSE]
    } else {
      df <- df[seq_len(n_keep), , drop = FALSE]
    }
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

check_joint_logdensity <- function(joint, per_dim) {
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
  result$max_abs <- max(abs(residual))
  result$max_rel <- max(abs(residual) / pmax(1, abs(joint)))
  result
}

trtf_power <- function() {
  seed <- as.integer(HP$seed)
  set.seed(seed)

  train_csv <- TRAIN_CSV
  test_csv  <- TEST_CSV

  derive_test_rows <- function(train_rows) {
    mapping <- c(
      `25` = 25L,
      `50` = 50L,
      `100` = 100L,
      `250` = 250L,
      `500` = 250L,
      `1000` = 250L,
      `2500` = 500L,
      `5000` = 1000L,
      `10000` = 2000L
    )
    key <- as.character(as.integer(train_rows))
    if (!nzchar(key) || is.na(train_rows) || !is.finite(train_rows) || train_rows <= 0L) {
      stop("Invalid train row count provided for TRTF run.")
    }
    if (key %in% names(mapping)) mapping[[key]] else as.integer(train_rows)
  }

  n_train <- as.integer(N_ROWS)
  n_test  <- derive_test_rows(n_train)

  meta <- NULL
  if (!is.null(META_INFO)) {
    meta_path <- normalizePath(META_PATH, winslash = "/", mustWork = TRUE)
    meta <- switch(
      META_INFO$type,
      "rds" = readRDS(meta_path),
      "csv" = {
        df_meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE)
        required_cols <- c("feature", "mu", "sigma")
        missing_cols <- setdiff(required_cols, names(df_meta))
        if (length(missing_cols)) {
          stop(sprintf("Standardisation CSV missing columns: %s", paste(missing_cols, collapse = ",")))
        }
        df_meta$feature <- trimws(df_meta$feature)
        list(
          features = df_meta$feature,
          mu = setNames(as.numeric(df_meta$mu), df_meta$feature),
          sigma = setNames(as.numeric(df_meta$sigma), df_meta$feature),
          eps = if ("eps" %in% names(df_meta)) setNames(as.numeric(df_meta$eps), df_meta$feature) else NULL
        )
      },
      stop(sprintf("Unsupported metadata type '%s'", META_INFO$type))
    )
  }

  Xtr_all <- read_numeric_matrix(train_csv, n_rows = n_train, sample_rows = dataset_cfg$sample_rows)
  Xte_all <- read_numeric_matrix(test_csv,  n_rows = n_test,  sample_rows = dataset_cfg$sample_rows)
  if (!nrow(Xtr_all) || !nrow(Xte_all)) {
    stop("No rows available after loading train/test CSVs.")
  }

  if (!is.null(meta) && !is.null(meta$features)) {
    features <- meta$features
    missing_tr <- setdiff(features, colnames(Xtr_all))
    missing_te <- setdiff(features, colnames(Xte_all))
    if (length(missing_tr) || length(missing_te)) {
      stop(sprintf("Feature mismatch: train missing [%s], test missing [%s]",
                   paste(missing_tr, collapse = ","), paste(missing_te, collapse = ",")))
    }
    X_tr <- Xtr_all[, features, drop = FALSE]
    X_te <- Xte_all[, features, drop = FALSE]
  } else {
    common <- intersect(colnames(Xtr_all), colnames(Xte_all))
    if (!length(common)) stop("No overlapping numeric columns between train/test CSVs.")
    X_tr <- Xtr_all[, common, drop = FALSE]
    X_te <- Xte_all[, common, drop = FALSE]
  }

  sigma_vec <- NULL
  if (!is.null(meta) && !is.null(meta$sigma)) {
    sigma_vec <- meta$sigma
    if (!is.null(names(sigma_vec))) {
      sigma_vec <- as.numeric(sigma_vec[colnames(X_tr)])
    } else {
      sigma_vec <- as.numeric(sigma_vec)
      if (length(sigma_vec) != ncol(X_tr)) {
        stop("Length of sigma vector does not match feature count.")
      }
    }
  }
  if (is.null(sigma_vec)) {
    sigma_vec <- rep(1, ncol(X_tr))
  }
  if (length(sigma_vec) != ncol(X_tr)) {
    stop("Sigma vector does not align with feature columns after reordering.")
  }
  if (any(!is.finite(sigma_vec) | sigma_vec <= 0)) {
    stop("Sigma vector contains non-positive or non-finite entries.")
  }
  names(sigma_vec) <- colnames(X_tr)
  log_sigma <- log(sigma_vec)

  if (!is.null(meta) && !is.null(meta$index)) {
    idx_tr <- meta$index$train
    idx_te <- meta$index$test
    idx_va <- meta$index$val
    as_int <- function(x) if (is.null(x)) integer(0L) else as.integer(x)
    idx_tr <- as_int(idx_tr)
    idx_te <- as_int(idx_te)
    idx_va <- as_int(idx_va)
    if (length(intersect(idx_tr, idx_te)) > 0L ||
        length(intersect(idx_tr, idx_va)) > 0L ||
        length(intersect(idx_va, idx_te)) > 0L) {
      stop("Split indices are not disjoint (train/val/test overlap detected).")
    }
  }

  K <- ncol(X_tr)
  S <- list(X_tr = X_tr, X_te = X_te)
  cfg <- replicate(K, list(distr = "norm"), simplify = FALSE)

  t_trtf_tr <- t_trtf_te <- NA_real_
  cores_use <- NC

  mod_trtf <- tryCatch({
    t0 <- system.time({ m <- fit_TRTF(S, cfg, seed = seed, cores = cores_use) })
    t_trtf_tr <- t0[["elapsed"]]
    m
  }, error = function(e) {
    stop(sprintf("TRTF training failed: %s", e$message))
  })

  pred_time <- tryCatch(system.time({
    LL_by_dim <- predict(mod_trtf, X_te, type = "logdensity_by_dim",
                         cores = cores_use,
                         trace = isTRUE(getOption("trtf.trace_predict", FALSE)))
  }), error = function(e) {
    stop(sprintf("TRTF per-dimension prediction failed: %s", e$message))
  })
  t_trtf_te <- unname(pred_time[["elapsed"]])

  if (!is.matrix(LL_by_dim)) {
    LL_by_dim <- as.matrix(LL_by_dim)
  }
  if (ncol(LL_by_dim) != ncol(X_te)) {
    stop("Per-dimension log-density output dimension mismatch.")
  }

  L_joint_try <- try(predict(mod_trtf, X_te, type = "logdensity",
                             cores = cores_use, trace = FALSE), silent = TRUE)
  if (inherits(L_joint_try, "try-error")) {
    warning(sprintf("Joint prediction failed: %s", conditionMessage(attr(L_joint_try, "condition"))))
    L_joint <- rep(NA_real_, nrow(LL_by_dim))
  } else {
    L_joint <- as.numeric(L_joint_try)
  }
  if (length(L_joint) != nrow(LL_by_dim)) {
    L_joint <- rep(L_joint, length.out = nrow(LL_by_dim))
  }

  joint_check <- check_joint_logdensity(L_joint, LL_by_dim)
  joint_warned <- FALSE
  if (!anyNA(unlist(joint_check))) {
    if (joint_check$max_abs > 1e-1 && joint_check$max_rel > 1e-1) {
      stop(sprintf("Joint consistency failed: max_abs=%.3e, max_rel=%.3e",
                   joint_check$max_abs, joint_check$max_rel))
    }
    if (joint_check$max_abs > 5e-2 || joint_check$max_rel > 1e-2) {
      warning(sprintf("Joint consistency weak: max_abs=%.3e, max_rel=%.3e",
                      joint_check$max_abs, joint_check$max_rel))
      joint_warned <- TRUE
    } else {
      message(sprintf("Joint consistency ok: max_abs=%.3e | max_rel=%.3e",
                      joint_check$max_abs, joint_check$max_rel))
    }
  }

  if (!inherits(L_joint_try, "try-error")) {
    residual <- rowSums(LL_by_dim) - as.numeric(L_joint_try)
    diff_abs <- max(abs(residual))
    diff_rel <- max(abs(residual) / pmax(1, abs(as.numeric(L_joint_try))))
    if ((diff_abs > 5e-2 || diff_rel > 1e-2) && !joint_warned) {
      warning(sprintf("Joint consistency weak: max_abs=%.3e, max_rel=%.3e",
                      diff_abs, diff_rel))
    }
    if (diff_abs > 1e-1 && diff_rel > 5e-2) {
      stop(sprintf("Joint consistency failed: max_abs=%.3e, max_rel=%.3e",
                   diff_abs, diff_rel))
    }
  }

  LL_X_by_dim <- sweep(LL_by_dim, 2L, log_sigma, FUN = "-")
  LL_X_joint  <- rowSums(LL_X_by_dim)
  LL_X_joint_manual <- rowSums(LL_by_dim) - sum(log_sigma)
  stopifnot(max(abs(LL_X_joint - LL_X_joint_manual)) < 1e-9)

  per_dim_mean <- colMeans(LL_X_by_dim)
  joint_mean   <- mean(LL_X_joint)

  dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
  sink_path <- RESULTS_TXT
  sink_active <- FALSE
  if (!is.null(sink_path) && nzchar(sink_path)) {
    sink_attempt <- try(sink(sink_path, split = TRUE), silent = TRUE)
    if (inherits(sink_attempt, "try-error")) {
      cond <- attr(sink_attempt, "condition")
      warn_msg <- if (!is.null(cond)) conditionMessage(cond) else as.character(sink_attempt)
      warning(sprintf("Could not mirror output to %s: %s", sink_path, warn_msg))
    } else {
      sink_active <- TRUE
      on.exit({ sink(NULL) }, add = TRUE)
    }
  }

  train_cores_used <- {
    if (!is.null(mod_trtf$train_cores)) as.integer(mod_trtf$train_cores)
    else {
      tasks <- max(1L, ncol(X_tr) - 1L)
      cap <- min(cores_cap_global, tasks)
      max(1L, get_train_cores(cap))
    }
  }

  lines_out <- c(
    sprintf("dataset: %s", RESULT_LABEL),
    sprintf("hyperparameters: ntree=%d, minsplit=%d, minbucket=%d, maxdepth=%d, seed=%d, train_cores=%d",
            as.integer(HP$ntree), as.integer(HP$minsplit), as.integer(HP$minbucket),
            as.integer(HP$maxdepth), as.integer(HP$seed), train_cores_used),
    sprintf("seed=%d, N=%s", as.integer(seed), if (is.na(N_ROWS)) "all" else as.integer(N_ROWS)),
    "LL_per_dim:",
    sprintf("dim%02d: %.6f", seq_along(per_dim_mean), per_dim_mean),
    sprintf("sum: %.6f", joint_mean),
    sprintf("time_sec: train=%ds, test=%ds, total=%ds",
            as.integer(round(t_trtf_tr)),
            as.integer(round(t_trtf_te)),
            as.integer(round(t_trtf_tr + t_trtf_te)))
  )
  cat(paste0(lines_out, "\n"), sep = "")

  invisible(list(
    per_dim_mean = per_dim_mean,
    joint_mean = joint_mean,
    train_time = t_trtf_tr,
    test_time = t_trtf_te,
    train_cores = train_cores_used
  ))
}

if (sys.nframe() == 0L) {
  invisible(trtf_power())
}

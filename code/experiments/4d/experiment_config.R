if (!exists("locate_repo_loader", inherits = TRUE)) {
  locate_repo_loader <- function() {
    detect_script_path <- function() {
      frames <- sys.frames()
      for (i in rev(seq_along(frames))) {
        fi <- frames[[i]]
        if (!is.null(fi$ofile)) {
          path <- tryCatch(normalizePath(fi$ofile, winslash = "/", mustWork = TRUE),
                          error = function(e) NA_character_)
          if (!is.na(path) && nzchar(path)) return(path)
        }
      }
      args <- commandArgs(trailingOnly = FALSE)
      file_arg <- args[grepl("^--file=", args)]
      if (length(file_arg)) {
        cand <- sub("^--file=", "", file_arg[1])
        path <- tryCatch(normalizePath(cand, winslash = "/", mustWork = TRUE),
                         error = function(e) NA_character_)
        if (!is.na(path) && nzchar(path)) return(path)
      }
      NA_character_
    }

    start_dirs <- character()
    script_path <- detect_script_path()
    if (!is.na(script_path) && nzchar(script_path)) {
      start_dirs <- c(start_dirs, dirname(script_path))
    }
    wd <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = FALSE),
                   error = function(e) getwd())
    start_dirs <- unique(c(start_dirs, wd))
    checked <- character()
    for (start in start_dirs) {
      cur <- start
      repeat {
        cur <- tryCatch(normalizePath(cur, winslash = "/", mustWork = FALSE),
                        error = function(e) cur)
        if (!nzchar(cur) || cur %in% checked) break
        checked <- c(checked, cur)
        cand1 <- file.path(cur, "R", "loader.R")
        if (file.exists(cand1)) {
          return(normalizePath(cand1, winslash = "/", mustWork = TRUE))
        }
        cand2 <- file.path(cur, "code", "R", "loader.R")
        if (file.exists(cand2)) {
          return(normalizePath(cand2, winslash = "/", mustWork = TRUE))
        }
        parent <- dirname(cur)
        if (identical(parent, cur)) break
        cur <- parent
      }
    }
    stop("Could not locate loader.R")
  }
}

loader_path <- locate_repo_loader()
if (!exists("initialize_repo")) {
  source(loader_path, chdir = FALSE)
}
root_path <- initialize_repo()

if (!exists("config", inherits = FALSE)) {
  stop("`config` must be defined before sourcing main_config.R")
}
if (!exists("n", inherits = FALSE)) {
  if (exists("N", inherits = FALSE)) {
    n <- N
  } else {
    stop("`n` (or `N`) must be defined before sourcing main_config.R")
  }
}

# Order heuristics now included in ttm_core.R

# RStudio-friendly power-run config (set enabled <- TRUE and paths, then source)
POWER_CONFIG <- list(
  enabled = FALSE,
  train_csv = file.path(root_path, "experiments", "normalizing flow", "power", "power_train_head_1pct.csv"),
  test_csv  = file.path(root_path, "experiments", "normalizing flow", "power", "power_test_head_1pct.csv"),
  n = 1000L,
  models = c("trtf", "ttm", "ttm_sep", "copula_np")
)

run_power_inline <- function(train_csv, test_csv, n_train = 1000L, n_test = 1000L, models = c("trtf","ttm","ttm_sep","copula_np")) {
  read_num <- function(path, n) {
    df <- utils::read.csv(path, header = TRUE, check.names = TRUE)
    num <- df[vapply(df, is.numeric, logical(1))]
    if (n > 0 && n < nrow(num)) num <- utils::head(num, n)
    as.matrix(num)
  }
  Xtr_all <- read_num(train_csv, n_train)
  Xte_all <- read_num(test_csv,  n_test)
  common <- intersect(colnames(Xtr_all), colnames(Xte_all))
  stopifnot(length(common) > 0)
  X_tr <- Xtr_all[, common, drop = FALSE]
  X_te <- Xte_all[, common, drop = FALSE]
  S <- list(X_tr = X_tr, X_te = X_te)
  mods <- list()
  if ("trtf" %in% models) mods$trtf <- tryCatch(fit_TRTF(S, replicate(ncol(X_tr), list(distr = "norm"), simplify = FALSE), seed = 42), error = function(e) NULL)
  if ("ttm" %in% models) mods$ttm <- fit_ttm(S, algo = "marginal", seed = 42)$S
  if ("ttm_sep" %in% models) mods$ttm_sep <- fit_ttm(S, algo = "separable", seed = 42)$S
  if ("copula_np" %in% models && exists("fit_copula_np")) mods$copula_np <- fit_copula_np(S, seed = 42)
  # Build simple table
  N <- nrow(X_te); K <- ncol(X_te)
  stderr <- function(x) stats::sd(x) / sqrt(length(x))
  fmt <- function(m, se) sprintf("%.2f ± %.2f", round(m, 2), round(2 * se, 2))
  tab <- data.frame(dim = as.character(seq_len(K)), distribution = rep("data", K), stringsAsFactors = FALSE)
  if (!is.null(mods$trtf)) { ll <- tryCatch(-predict(mods$trtf, X_te, type = "logdensity_by_dim"), error = function(e) matrix(NA_real_, N, K)); tab[["Random Forest"]] <- fmt(colMeans(ll), apply(ll, 2, stderr)) }
  if (!is.null(mods$ttm))  { ll <- -predict(mods$ttm, X_te, type = "logdensity_by_dim"); tab[["Marginal Map"]] <- fmt(colMeans(ll), apply(ll, 2, stderr)) }
  if (!is.null(mods$ttm_sep)) { ll <- -predict(mods$ttm_sep, X_te, type = "logdensity_by_dim"); tab[["Separable Map"]] <- fmt(colMeans(ll), apply(ll, 2, stderr)) }
  if (!is.null(mods$copula_np)) { ll <- tryCatch(-predict(mods$copula_np, X_te, type = "logdensity_by_dim"), error = function(e) matrix(NA_real_, N, K)); tab[["Copula NP"]] <- fmt(colMeans(ll), apply(ll, 2, stderr)) }
  sum_row <- data.frame(dim = "k", distribution = "SUM", stringsAsFactors = FALSE)
  if (!is.null(mods$trtf))     sum_row[["Random Forest"]] <- fmt(sum(colMeans(tryCatch(-predict(mods$trtf, X_te, type = "logdensity_by_dim"), error=function(e) matrix(NA_real_,N,K)))), NA_real_)
  if (!is.null(mods$ttm))      sum_row[["Marginal Map"]]  <- fmt(sum(colMeans(-predict(mods$ttm, X_te, type = "logdensity_by_dim"))), NA_real_)
  if (!is.null(mods$ttm_sep))  sum_row[["Separable Map"]] <- fmt(sum(colMeans(-predict(mods$ttm_sep, X_te, type = "logdensity_by_dim"))), NA_real_)
  if (!is.null(mods$copula_np))sum_row[["Copula NP"]]     <- fmt(sum(colMeans(tryCatch(-predict(mods$copula_np, X_te, type = "logdensity_by_dim"), error=function(e) matrix(NA_real_,N,K)))), NA_real_)
  tab <- rbind(tab, sum_row)
  cat(sprintf("Train CSV: %s\nTest CSV: %s\n", train_csv, test_csv))
  cat(sprintf("n_train=%d, n_test=%d, K=%d\n", n_train, n_test, ncol(X_tr)))
  print(tab)
  invisible(tab)
}


perm <- c(1, 2, 3, 4)

#' @export
main <- function() {
  # If POWER_CONFIG enabled, run power inline and return
  if (is.list(POWER_CONFIG) && isTRUE(POWER_CONFIG$enabled)) {
    return(invisible(run_power_inline(POWER_CONFIG$train_csv, POWER_CONFIG$test_csv, as.integer(POWER_CONFIG$n), as.integer(POWER_CONFIG$n), POWER_CONFIG$models)))
  }
  seed <- as.integer(Sys.getenv("SEED", 42))
  set.seed(seed)
  prep <- prepare_data(n, config, seed = seed)
  S0 <- prep$S
  # Use explicit permutation defined at top-level 'perm'
  if (!is.numeric(perm) || length(perm) != ncol(S0$X_tr)) {
    perm <<- seq_len(ncol(S0$X_tr))
  }
  message(sprintf("[ORDER] Using explicit permutation: %s", paste(perm, collapse = ",")))
  S <- list(
    X_tr  = S0$X_tr[, perm, drop = FALSE],
    X_te  = S0$X_te[, perm, drop = FALSE]
  )
  cfg <- config[perm]
  
  t_true_tr  <- system.time(mod_true      <- fit_TRUE(S, cfg))[['elapsed']]
  t_joint_tr <- 0
  mod_true_joint <- tryCatch({
    fit_TRUE_JOINT(S, cfg)
  }, error = function(e) {
    message("[WARN] Skipping True (Joint) for this permutation: ", e$message)
    NULL
  })
  t_trtf_tr  <- system.time(mod_trtf      <- fit_TRTF(S, cfg, seed = seed))[['elapsed']]
  # Copula NP baseline (robust: 2D-labeled copula or KDE product fallback)
  t_cop_tr   <- tryCatch(system.time(mod_cop      <- fit_copula_np(S, seed = seed))[["elapsed"]], error = function(e) NA_real_)
  # New TTM fits (maps-from-samples, forward-KL)
  mod_ttm      <- fit_ttm(S, algo = "marginal",  seed = seed);  t_ttm_tr <- mod_ttm$time_train
  mod_ttm_sep  <- fit_ttm(S, algo = "separable", seed = seed);  t_sep_tr <- mod_ttm_sep$time_train
  # Cross-term removed

  t_true_te  <- system.time(logL_TRUE(mod_true, S$X_te))[['elapsed']]
  t_joint_te <- tryCatch({
    system.time(true_joint_logdensity_by_dim(config, S0$X_te))[['elapsed']]
  }, error = function(e) NA_real_)
  t_trtf_te  <- tryCatch(system.time(predict(mod_trtf, S$X_te, type = "logdensity_by_dim"))[['elapsed']], error = function(e) NA_real_)
  t_cop_te   <- tryCatch(system.time(predict(mod_cop, S$X_te, type = "logdensity_by_dim"))[['elapsed']], error = function(e) NA_real_)
  t_ttm_te   <- system.time(predict_ttm(mod_ttm$S,       S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_sep_te   <- system.time(predict_ttm(mod_ttm_sep$S,   S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  

  mods <- list(
    true = mod_true,
    true_joint = mod_true_joint,
    trtf = tryCatch(mod_trtf, error = function(e) NULL),
    ttm  = mod_ttm,
    ttm_sep = mod_ttm_sep
  )
  tab <- calc_loglik_tables(mods, cfg, S$X_te, config_canonical = config, perm = perm)
  if ("train_test_policy" %in% names(tab)) {
    tab$train_test_policy <- NULL
  }
  fmt <- function(m, se) sprintf("%.2f ± %.2f", round(m, 2), round(2 * se, 2))
  # Append Copula NP column using the unified predict API
  LD_cop <- tryCatch(predict(mod_cop, S$X_te, type = "logdensity_by_dim"), error = function(e) NULL)
  if (!is.null(LD_cop) && is.matrix(LD_cop) && all(dim(LD_cop) == dim(S$X_te))) {
    per_cop <- -colMeans(LD_cop)
    se_cop  <- apply(-LD_cop, 2, stderr)
    tab[["Copula NP"]] <- c(
      fmt(per_cop, se_cop),
      fmt(sum(per_cop), stats::sd(rowSums(-LD_cop)) / sqrt(nrow(S$X_te)))
    )
  } else {
    tab[["Copula NP"]] <- c(rep("NA", ncol(S$X_te)), "NA")
  }
  # no composite column
  cat(sprintf("n=%d\n", n))
  print(tab)
  cat(sprintf("Permutation order %s (train/test only)\n", paste(perm, collapse = ",")))
  time_tab <- data.frame(
    model = c("True (marginal)", "True (Joint)", "Random Forest",
              "Copula NP", "Marginal Map", "Separable Map"),
    train_sec = c(t_true_tr, t_joint_tr, t_trtf_tr,
                  t_cop_tr, t_ttm_tr, t_sep_tr),
    test_sec = c(t_true_te, t_joint_te, t_trtf_te,
                 t_cop_te, t_ttm_te, t_sep_te),
    stringsAsFactors = FALSE
  )
  time_tab$total_sec <- with(time_tab, train_sec + test_sec)
  stopifnot(all.equal(time_tab$total_sec,
                      time_tab$train_sec + time_tab$test_sec))
  print(time_tab)
  timing_table <<- time_tab
  results_table <<- tab
  stopifnot(identical(results_table, tab))
  return(tab)
}

if (sys.nframe() == 0L) {
  invisible(main())
}

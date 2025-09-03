get_repo_root <- function() {
  p <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:10) {
    if (file.exists(file.path(p, "ALGORITHM_SPEC.md")) && file.exists(file.path(p, "00_globals.R"))) return(p)
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
.repo_root <- if (exists("root_path")) root_path else get_repo_root()
source(file.path(.repo_root, "00_globals.R"))
source(file.path(.repo_root, "01_data_generation.R"))
source(file.path(.repo_root, "02_split.R"))
source(file.path(.repo_root, "models/true_model.R"))
source(file.path(.repo_root, "models/trtf_model.R"))
source(file.path(.repo_root, "models/ttm/ttm_bases.R"))
source(file.path(.repo_root, "models/ttm/ttm_core.R"))
source(file.path(.repo_root, "models/ttm/ttm_marginal.R"))
source(file.path(.repo_root, "models/ttm/ttm_separable.R"))
source(file.path(.repo_root, "models/ttm/ttm_crossterm.R"))
source(file.path(.repo_root, "models/true_joint_model.R"))
source(file.path(.repo_root, "04_evaluation.R"))
## Optional: order heuristics
if (file.exists(file.path(.repo_root, "models", "ttm", "order_heuristics.R"))) {
  source(file.path(.repo_root, "models", "ttm", "order_heuristics.R"))
}
## Optional nonparametric copula baseline
if (file.exists(file.path(.repo_root, "models", "copula_np.R"))) {
  source(file.path(.repo_root, "models", "copula_np.R"))
}


perm <- c(1, 2, 3, 4)
n <- 50
# Optional: allow overriding n via environment variable N_OVERRIDE
nv <- suppressWarnings(as.integer(Sys.getenv("N_OVERRIDE", "")))
if (is.finite(nv) && !is.na(nv) && nv > 0L) n <- nv
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "beta",
       parm = function(d) list(shape1 = softplus(d$X2),
                               shape2 = softplus(d$X1))),
  list(distr = "gamma",
       parm = function(d) list(shape = softplus(d$X3),
                               scale = softplus(d$X2)))
)

#' @export
main <- function() {
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
  # Cross-term with moderate expressiveness and accuracy
  # CPU-freundliche Defaults: Auto-λ, moderate Basis, stärkere Sparsity
  df_t_opt <- as.integer(max(16L, floor(n / 10)))
  options(
    cross.lambda_non = NULL, cross.lambda_mon = NULL,
    cross.df_t = 12L,
    cross.sparsity_tau = 0.10
  )
  mod_ttm_cross <- fit_ttm(S, algo = "crossterm", seed = seed,
                           deg_g = 2L, df_t = 6L, Q = 12L,
                           lambda = NA_real_, Hmax = 4L, maxit = 180L,
                           order_search = list(
                             strategy   = "greedy",
                             start      = "chol_pivot",
                             val_frac   = 0.2,
                             max_orders = 8L,
                             maxit_fast = 35L,
                             Q_fast     = 8L,
                             eps        = 1e-3
                           ))
  t_ct_tr  <- mod_ttm_cross$time_train

  t_true_te  <- system.time(logL_TRUE(mod_true, S$X_te))[['elapsed']]
  t_joint_te <- tryCatch({
    system.time(true_joint_logdensity_by_dim(config, S0$X_te))[['elapsed']]
  }, error = function(e) NA_real_)
  t_trtf_te  <- tryCatch(system.time(predict(mod_trtf, S$X_te, type = "logdensity_by_dim"))[['elapsed']], error = function(e) NA_real_)
  t_cop_te   <- tryCatch(system.time(predict(mod_cop, S$X_te, type = "logdensity_by_dim"))[['elapsed']], error = function(e) NA_real_)
  t_ttm_te   <- system.time(predict_ttm(mod_ttm$S,       S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_sep_te   <- system.time(predict_ttm(mod_ttm_sep$S,   S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_ct_te    <- system.time(predict_ttm_crossterm(mod_ttm_cross$S, S$X_te, type = "logdensity_by_dim"))[['elapsed']]

  mods <- list(
    true = mod_true,
    true_joint = mod_true_joint,
    trtf = tryCatch(mod_trtf, error = function(e) NULL),
    ttm  = mod_ttm,
    ttm_sep = mod_ttm_sep,
    ttm_cross = mod_ttm_cross
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
  cat(sprintf("n=%d\n", n))
  print(tab)
  cat(sprintf("Permutation order %s (train/test only)\n", paste(perm, collapse = ",")))
  time_tab <- data.frame(
    model = c("True (marginal)", "True (Joint)", "Random Forest",
              "Copula NP", "Marginal Map", "Separable Map", "Cross-term Map"),
    train_sec = c(t_true_tr, t_joint_tr, t_trtf_tr,
                  t_cop_tr, t_ttm_tr, t_sep_tr, t_ct_tr),
    test_sec = c(t_true_te, t_joint_te, t_trtf_te,
                 t_cop_te, t_ttm_te, t_sep_te, t_ct_te),
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

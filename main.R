source("00_globals.R")
if (!exists("%||%")) "%||%" <- function(a, b) if (is.null(a)) b else a
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ttm_marginal.R")
source("models/ttm_separable.R")
source("models/ttm_cross_term.R")
source("models/true_joint_model.R")
source("04_evaluation.R")

perm <- c(1, 2, 3, 4)
n <- 50
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

#' Starte die komplette Analyse
#' 
#' Es werden zwei Log-Likelihood-Tabellen ausgegeben: einmal für die
#' normale Reihenfolge der Dimensionen und einmal für eine vorgegebene
#' Permutation.
#' 
#' @export
main <- function() {
  seed <- as.integer(Sys.getenv("SEED", 42))
  set.seed(seed)
  prep <- prepare_data(n, config, seed = seed)
  S0 <- prep$S
  S <- list(
    X_tr  = S0$X_tr[, perm, drop = FALSE],
    X_te  = S0$X_te[, perm, drop = FALSE]
  )
  cfg <- config[perm]

  t_true_tr  <- system.time(mod_true      <- fit_TRUE(S, cfg))[['elapsed']]
  t_joint_tr <- 0
  mod_true_joint <- fit_TRUE_JOINT(S, cfg)
  t_trtf_tr  <- system.time(mod_trtf      <- fit_TRTF(S, cfg, seed = seed))[['elapsed']]
  mod_ttm     <- trainMarginalMap(S, seed = seed);   t_ttm_tr <- mod_ttm$time_train
  mod_ttm_sep <- trainSeparableMap(S, seed = seed);  t_sep_tr <- mod_ttm_sep$time_train
  # Choose reverse-KL to NF if enabled, else forward-KL
  use_rev_nf <- Sys.getenv('MDE_CTM_USE_REV_NF','0') %in% c('1','true','TRUE')
  if (use_rev_nf) {
    mod_ttm_cross <- trainCrossTermMapReverseNF(S, degree_g = 2, degree_t = 2, seed = seed)
    t_ct_tr <- NA_real_
  } else {
    mod_ttm_cross <- trainCrossTermMap(S, degree_g = 3, seed = seed, warmstart_from_separable = TRUE); t_ct_tr <- mod_ttm_cross$time_train
  }
  # Log ridge configuration for TTM-Cross
  if (!is.null(mod_ttm_cross$S$meta$ridge)) {
    rr <- mod_ttm_cross$S$meta$ridge
    cat(sprintf("[RIDGE] lambda_non=%.3g, lambda_mon=%.3g\n", rr[["lambda_non"]], rr[["lambda_mon"]]))
  } else if (!is.null(mod_ttm_cross$S$meta$lambda_non) && !is.null(mod_ttm_cross$S$meta$lambda_mon)) {
    cat(sprintf("[RIDGE] lambda_non=%.3g, lambda_mon=%.3g\n", mod_ttm_cross$S$meta$lambda_non, mod_ttm_cross$S$meta$lambda_mon))
  }

  t_true_te  <- system.time(logL_TRUE(mod_true, S$X_te))[['elapsed']]
  t_joint_te <- system.time(true_joint_logdensity_by_dim(cfg, S$X_te))[['elapsed']]
  t_trtf_te  <- system.time(predict(mod_trtf, S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_ttm_te   <- mod_ttm$time_pred
  t_sep_te   <- mod_ttm_sep$time_pred
  t_ct_te    <- mod_ttm_cross$time_pred

  mods <- list(
    true = mod_true,
    true_joint = mod_true_joint,
    trtf = mod_trtf,
    ttm  = mod_ttm,
    ttm_sep = mod_ttm_sep,
    ttm_cross = mod_ttm_cross
  )
  tab <- calc_loglik_tables(mods, cfg, S$X_te)
  cat(sprintf("n=%d\n", n))
  print(tab)
  cat(sprintf("Permutation order %s (train/test only)\n", paste(perm, collapse = ",")))
  time_tab <- data.frame(
    model = c("True (marginal)", "True (Joint)", "Random Forest",
              "Marginal Map", "Separable Map", "Cross-term Map"),
    train_sec = c(t_true_tr, t_joint_tr, t_trtf_tr,
                  t_ttm_tr, t_sep_tr, t_ct_tr),
    test_sec = c(t_true_te, t_joint_te, t_trtf_te,
                 t_ttm_te, t_sep_te, t_ct_te),
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
  main()
}

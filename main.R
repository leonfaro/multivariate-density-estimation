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
source("replicate_code.R")

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
  set.seed(42)
  prep <- prepare_data(n, config, seed = 42)
  S0 <- prep$S
  S <- list(
    X_tr  = S0$X_tr[, perm, drop = FALSE],
    X_val = S0$X_val[, perm, drop = FALSE],
    X_te  = S0$X_te[, perm, drop = FALSE]
  )
  cfg <- config[perm]

  t_true_tr  <- system.time(mod_true      <- fit_TRUE(S, cfg))[['elapsed']]
  t_joint_tr <- 0
  mod_true_joint <- fit_TRUE_JOINT(S, cfg)
  t_trtf_tr  <- system.time(mod_trtf      <- fit_TRTF(S, cfg, seed = 42))[['elapsed']]
  mod_ttm     <- trainMarginalMap(S);  t_ttm_tr <- mod_ttm$time_train
  mod_ttm_sep <- trainSeparableMap(S); t_sep_tr <- mod_ttm_sep$time_train
  mod_ttm_cross <- trainCrossTermMap(S); t_ct_tr <- mod_ttm_cross$time_train

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
  cat(sprintf("Permutation order %s\n", paste(perm, collapse = ",")))
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
  invisible(tab)
}

if (sys.nframe() == 0L) {
  main()
  replicate_code_scripts("main.R", "replicated_code.txt", env = globalenv())
}


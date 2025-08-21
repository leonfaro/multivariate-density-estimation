source("00_globals.R")
if (!exists("%||%")) "%||%" <- function(a, b) if (is.null(a)) b else a
source("01_data_generation.R")
source("02_split.R")
source("scripts/halfmoon_data.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ttm_marginal.R")
source("models/ttm_separable.R")
source("models/ttm_cross_term.R")
source("models/true_joint_model.R")
source("04_evaluation.R")
source("replicate_code.R")

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
  dataset <- Sys.getenv("DATASET", "config4d")
  seed <- as.integer(Sys.getenv("SEED", 42))
  set.seed(seed)

  if (dataset == "halfmoon2d") {
    ntr <- pmin(as.integer(Sys.getenv("N_TRAIN", 250)), 250)
    nte <- pmin(as.integer(Sys.getenv("N_TEST", 250)), 250)
    noise <- as.numeric(Sys.getenv("NOISE", 0.15))
    S <- make_halfmoon_splits(ntr, nte, noise, seed)
    S$meta$dataset <- dataset
    dir.create("results", showWarnings = FALSE)
    saveRDS(S, sprintf("results/splits_%s_seed%03d.rds", dataset, seed))
    cat(sprintf("[DATASET %s] K=%d | n_tr=%d (val=%d) | n_te=%d | noise=%.3f | seed=%d\n",
                dataset, ncol(S$X_tr), nrow(S$X_tr), nrow(S$X_val),
                nrow(S$X_te), noise, seed))
    set.seed(seed)
    mods <- list()
    eval_halfmoon(mods, S, NULL)
    invisible(NULL)
  } else {
    prep <- prepare_data(n, config, seed = seed)
    mods <- list(
      true = fit_TRUE(prep$S, config),
      true_joint = fit_TRUE_JOINT(prep$S, config),
      trtf = fit_TRTF(prep$S, config, seed = seed),
      ttm  = trainMarginalMap(prep$S),
      ttm_sep = trainSeparableMap(prep$S),
      ttm_cross = trainCrossTermMap(prep$S)
    )
    tab <- calc_loglik_tables(mods, config, prep$S$X_te)
    print(tab)
    results_table <<- tab
    invisible(tab)
  }
}

if (sys.nframe() == 0L) {
  main()
  replicate_code_scripts("main.R", "replicated_code.txt", env = globalenv())
}


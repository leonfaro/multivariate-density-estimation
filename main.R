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
  mods <- list(
    true = fit_TRUE(S, cfg),
    true_joint = fit_TRUE_JOINT(S, cfg),
    trtf = fit_TRTF(S, cfg, seed = 42),
    ttm  = trainMarginalMap(S),
    ttm_sep = trainSeparableMap(S),
    ttm_cross = trainCrossTermMap(S)
  )
  tab <- calc_loglik_tables(mods, cfg, S$X_te)
  cat(sprintf("n=%d\n", n))
  print(tab)
  results_table <<- tab
  invisible(tab)
}

if (sys.nframe() == 0L) {
  main()
  replicate_code_scripts("main.R", "replicated_code.txt", env = globalenv())
}


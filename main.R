source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ttm_marginal.R")
source("models/ttm_separable.R")
source("models/ttm_cross_term.R")
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
  set.seed(42)
  prep <- prepare_data(n, config, seed = 42)
  mods <- list(
    true = fit_TRUE(prep$S, config),
    trtf = fit_TRTF(prep$S, config, seed = 42),
    ttm  = trainMarginalMap(prep$S),
    ttm_sep = trainSeparableMap(prep$S),
    ttm_cross = trainCrossTermMap(prep$S)
  )
  tab <- calc_loglik_tables(mods, config, prep$S$X_te)
  print(tab)
  results_table <<- tab
  invisible(tab)
}

if (sys.nframe() == 0L) {
  main()
  replicate_code_scripts("main.R", "replicated_code.txt", env = globalenv())
}


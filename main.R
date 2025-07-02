source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ks_model.R")
source("04_evaluation.R")

n <- 50
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "beta", parm = function(d) list(shape1 = softplus(d$X2), shape2 = softplus(d$X1))),
  list(distr = "gamma", parm = function(d) list(shape = softplus(d$X3), scale = softplus(d$X2)))
)

#' Starte die komplette Analyse
#' 
#' Es werden zwei Log-Likelihood-Tabellen ausgegeben: einmal für die
#' normale Reihenfolge der Dimensionen und einmal für eine vorgegebene
#' Permutation.
#' 
#' @export
main <- function() {
  prep <- prepare_data(n, config, seed = 42)
  mods <- list(
    true = fit_TRUE(prep$S, config),
    trtf = fit_TRTF(prep$S, config),
    ks   = fit_KS(prep$S, config)
  )
  tab <- calc_loglik_tables(mods, config, prep$S$X_te)
  print(tab)
  invisible(tab)
}

if (sys.nframe() == 0L) {
  main()
}

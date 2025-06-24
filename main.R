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
perm <- c(3, 4, 1, 2)

#' Starte die komplette Analyse
#'
#' Es werden zwei Log-Likelihood-Tabellen ausgegeben: einmal für die
#' normale Reihenfolge der Dimensionen und einmal für eine vorgegebene
#' Permutation.
#'
#' @export
main <- function() {
  prep <- prepare_data(n, config, perm)

  mods_norm <- fit_models(prep$S, config)
  tab_normal <- calc_loglik_tables(mods_norm, config)

  mods_perm <- fit_models(prep$S_perm, config)
  tab_perm <- calc_loglik_tables(mods_perm, config)

  print(tab_normal)
  print(tab_perm)
  invisible(list(normal = tab_normal, perm = tab_perm))
}

if (sys.nframe() == 0L) {
  main()
}

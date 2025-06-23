source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ks_model.R")
source("04_evaluation.R")
source("EDA.R")

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
#' Diese Wrapper-Funktion ruft `run_pipeline()` aus `EDA.R` auf und gibt
#' die erzeugte Ergebnistabelle zurück.
#'
#' @export
main <- function() {
  run_pipeline(n, config, perm)
}

if (sys.nframe() == 0L) {
  main()
}

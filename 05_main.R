#' Orchestrate full experiment workflow
#'
#' This script assembles all individual steps defined in the repository.
#' It relies on the helper functions from Scripts 1--6 as documented in
#' `roadmap.md`.
#'
#' The lines below allow overriding the global experiment size `N`
#' and the distribution configuration `config` without modifying
#' `00_globals.R`.

source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("04_evaluation.R")

N <- 500
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = d$X1)),
  list(distr = "beta", parm = function(d) list(shape1 = d$X2, shape2 = 1)),
  list(distr = "gamma", parm = function(d) list(shape = d$X3, scale = 1))
)

#' @export
main <- function() {
  G <- list(
    N = N,
    config = config,
    seed = 42,
    split_ratio = 0.70
  )

  X <- gen_samples(G)                             # Script 2
  S <- train_test_split(X, G$split_ratio, G$seed) # Script 3
  M_TRUE <- fit_TRUE(S$X_tr, S$X_te, G$config)    # Script 5
  models <- setNames(list(M_TRUE), "TRUE")
  results <- evaluate_all(S$X_te, models)         # Script 6
  baseline_ll <- logL_TRUE_dim(M_TRUE, S$X_te)
  tab <- data.frame(
    dim = as.character(seq_along(G$config)),
    distribution = sapply(G$config, `[[`, "distr"),
    logL_baseline = baseline_ll,
    stringsAsFactors = FALSE
  )

  sum_row <- setNames(vector("list", ncol(tab)), names(tab))
  for (nm in names(tab)) {
    if (nm == "dim") {
      sum_row[[nm]] <- "k"
    } else if (is.numeric(tab[[nm]])) {
      sum_row[[nm]] <- sum(tab[[nm]])
    } else {
      sum_row[[nm]] <- NA
    }
  }
  tab <- rbind(tab, as.data.frame(sum_row, stringsAsFactors = FALSE))
  print(tab)
  invisible(tab)
}

# einfache Möglichkeit, weitere Modelle anzuhängen:
# - Quellcode in models/<neues>.R mit
#     fit_<NAME>()  und  logL_<NAME>()
# - anschließend hier laden und der Liste `models` hinzufügen.

if (sys.nframe() == 0L) {
  main()
}

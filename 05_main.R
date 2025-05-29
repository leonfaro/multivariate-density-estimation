#' Orchestrate full experiment workflow
#'
#' This script assembles all individual steps defined in the repository.
#' It relies on the helper functions from Scripts 1--6 as documented in
#' `roadmap.md`.
#'
#' @export
main <- function() {
  G <- setup_global()                             # Script 1
  X <- gen_samples(G)                             # Script 2
  S <- train_test_split(X, G$split_ratio, G$seed) # Script 3
  M_TTM  <- fit_TTM(S$X_tr, S$X_te, G$H_grid)     # Script 4
  M_TRUE <- fit_TRUE(S$X_tr, S$X_te, G$config)    # Script 5
  models <- list(TTM = M_TTM, TRUE = M_TRUE)
  results <- evaluate_all(S$X_te, models)         # Script 6
  print(results)
  invisible(results)
}

# einfache Möglichkeit, weitere Modelle anzuhängen:
# - Quellcode in models/<neues>.R mit
#     fit_<NAME>()  und  logL_<NAME>()
# - anschließend hier laden und der Liste `models` hinzufügen.

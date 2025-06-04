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
source("models/trtf_model.R")
source("models/ks_model.R")
source("04_evaluation.R")

N <- 100 
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "beta", parm = function(d) list(shape1 = softplus(d$X2), shape2 = softplus(d$X1))),
  list(distr = "gamma", parm = function(d) list(shape = softplus(d$X3), scale = softplus(d$X2)))
)
perm <- c(3, 4, 1, 2)

#' @export
main <- function() {
  G <- list(
    N = N,
    config = config,
    seed = 42,
    split_ratio = 0.5
  )

  X <- gen_samples(G)                             # Script 2
  S <- train_test_split(X, G$split_ratio, G$seed) # Script 3

  ## Modelle in Originalreihenfolge
  M_TRUE <- fit_TRUE(S$X_tr, S$X_te, G$config)    # Script 5
  M_TRTF <- fit_TRTF(S$X_tr, S$X_te, G$config)    # Script models/trtf_model.R
  M_KS   <- fit_KS(S$X_tr, S$X_te, G$config)      # Script models/ks_model.R

  baseline_ll <- logL_TRUE_dim(M_TRUE, S$X_te)
  trtf_ll     <- logL_TRTF_dim(M_TRTF, S$X_te)
  ks_ll       <- logL_KS_dim(M_KS, S$X_te)

  tab_normal <- data.frame(
    dim = as.character(seq_along(G$config)),
    distribution = sapply(G$config, `[[`, "distr"),
    logL_baseline = baseline_ll,
    logL_trtf = trtf_ll,
    logL_ks = ks_ll,
    stringsAsFactors = FALSE
  )

  sum_row <- setNames(vector("list", ncol(tab_normal)), names(tab_normal))
  for (nm in names(tab_normal)) {
    if (nm == "dim") {
      sum_row[[nm]] <- "k"
    } else if (is.numeric(tab_normal[[nm]])) {
      sum_row[[nm]] <- sum(tab_normal[[nm]])
    } else {
      sum_row[[nm]] <- NA
    }
  }
  tab_normal <- rbind(tab_normal, as.data.frame(sum_row, stringsAsFactors = FALSE))

  ## Modelle mit Permutation
  X_tr_p <- S$X_tr[, perm, drop = FALSE]
  X_te_p <- S$X_te[, perm, drop = FALSE]
  colnames(X_tr_p) <- paste0("X", seq_len(ncol(X_tr_p)))
  colnames(X_te_p) <- paste0("X", seq_len(ncol(X_te_p)))

  M_TRUE_p <- fit_TRUE(X_tr_p, X_te_p, G$config)
  M_TRTF_p <- fit_TRTF(X_tr_p, X_te_p, G$config)
  M_KS_p   <- fit_KS(X_tr_p, X_te_p, G$config)

  baseline_ll_p <- logL_TRUE_dim(M_TRUE_p, X_te_p)
  trtf_ll_p     <- logL_TRTF_dim(M_TRTF_p, X_te_p)
  ks_ll_p       <- logL_KS_dim(M_KS_p, X_te_p)

  tab_perm <- data.frame(
    dim = as.character(seq_along(G$config)),
    distribution = sapply(G$config, `[[`, "distr"),
    logL_baseline = baseline_ll_p,
    logL_trtf = trtf_ll_p,
    logL_ks = ks_ll_p,
    stringsAsFactors = FALSE
  )

  sum_row <- setNames(vector("list", ncol(tab_perm)), names(tab_perm))
  for (nm in names(tab_perm)) {
    if (nm == "dim") {
      sum_row[[nm]] <- "k"
    } else if (is.numeric(tab_perm[[nm]])) {
      sum_row[[nm]] <- sum(tab_perm[[nm]])
    } else {
      sum_row[[nm]] <- NA
    }
  }
  tab_perm <- rbind(tab_perm, as.data.frame(sum_row, stringsAsFactors = FALSE))

  tab_normal_print <- tab_normal
  num_cols <- sapply(tab_normal_print, is.numeric)
  tab_normal_print[num_cols] <- lapply(tab_normal_print[num_cols], round, 3)

  tab_perm_print <- tab_perm
  num_cols <- sapply(tab_perm_print, is.numeric)
  tab_perm_print[num_cols] <- lapply(tab_perm_print[num_cols], round, 3)

  print(tab_normal_print)
  print(tab_perm_print)
  invisible(list(normal = tab_normal, permutation = tab_perm))
}

# einfache Möglichkeit, weitere Modelle anzuhängen:
# - Quellcode in models/<neues>.R mit
#     fit_<NAME>()  und  logL_<NAME>()
# - anschließend hier laden und der Liste `models` hinzufügen.

if (sys.nframe() == 0L) {
  main()
}


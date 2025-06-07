#' Evaluate models on a test set
#'
#' The function collects negative log-likelihoods for a list of fitted models and
#' returns them ordered by performance. Model identifiers are inferred from the
#' names of `model_list` and translated to their respective `logL_<ID>` function.
#'
#' @param X_te matrix of test observations
#' @param model_list named list of fitted model objects
#' @return data.frame with columns `model_id` and `neg_logL`
#' @export
evaluate_all <- function(X_te, model_list) {
  stopifnot(is.matrix(X_te), is.list(model_list))
  ids <- names(model_list)
  if (is.null(ids)) ids <- paste0("M", seq_along(model_list))

  P <- data.frame(model_id = ids, neg_logL = NA_real_,
                  stringsAsFactors = FALSE)

  for (i in seq_along(model_list)) {
    id <- ids[i]
    M <- model_list[[i]]
    fun_name <- paste0("logL_", id)
    if (!exists(fun_name, mode = "function")) {
      stop(sprintf("No log-likelihood function %s found", fun_name))
    }
    loss <- do.call(fun_name, list(M, X_te))
    if (!is.finite(loss)) stop("log-likelihood not finite")
    P$neg_logL[i] <- loss
  }

  P <- P[order(P$neg_logL), ]
  rownames(P) <- NULL
  P
}

#' Merge log-likelihood tables and add runtime information
#'
#' Combines results from the normal variable order and a permutation
#' of the variables.  Runtimes for the permutation order are computed
#' on the fly.  Returns a `kableExtra` object for direct inclusion in
#' LaTeX or HTML reports.
#'
#' @param tab_normal table created in `main()` using the original order
#' @param tab_perm   table created in `main()` on permuted variables
#' @param M_TRUE_p   TRUE model fitted on permuted data
#' @param M_TRTF_p   TRTF model fitted on permuted data
#' @param M_KS_p     KS model fitted on permuted data
#' @param M_TTM_p    TTM model fitted on permuted data
#' @param X_te_p     test matrix corresponding to the permutation
#' @param t_normal   named numeric vector with runtimes of the normal
#'                   order (names: true, trtf, ks, ttm)
#' @return `kableExtra` table with average log-likelihoods and runtimes
#' @export
combine_logL_tables <- function(tab_normal, tab_perm,
                                M_TRUE_p, M_TRTF_p, M_KS_p, M_TTM_p, X_te_p,
                                t_normal) {
  if (!requireNamespace("kableExtra", quietly = TRUE))
    install.packages("kableExtra", repos = "https://cloud.r-project.org")

  library(dplyr)
  library(tibble)
  library(kableExtra)

  t_true_p <- system.time(logL_TRUE_dim(M_TRUE_p, X_te_p))["elapsed"]
  t_trtf_p <- system.time(logL_TRTF_dim(M_TRTF_p, X_te_p))["elapsed"]
  t_ks_p   <- system.time(logL_KS_dim(M_KS_p,   X_te_p))["elapsed"]
  t_ttm_p  <- system.time(logL_TTM_dim(M_TTM_p, X_te_p))["elapsed"]

  clean_cols <- function(df, suffix) {
    df %>%
      rename(
        distr = distribution,
        true  = logL_baseline,
        trtf  = logL_trtf,
        ks    = logL_ks,
        ttm   = logL_ttm
      ) %>%
      rename_with(~ paste0(.x, "_", suffix),
                  c(true, trtf, ks, ttm))
  }

  tab_norm <- clean_cols(tab_normal, "norm")
  tab_perm <- clean_cols(tab_perm,   "perm")

  tab_all <- left_join(tab_norm, tab_perm, by = c("dim", "distr"))
  tab_all <- tab_all %>%
    mutate(across(where(is.numeric), replace_na, 0)) %>%
    mutate(distr = replace_na(distr, "SUM"))

  tab_all <- tab_all %>%
    add_row(
      dim        = "runtime",
      distr      = "",
      true_norm  = t_normal["true"] * 1000,
      true_perm  = t_true_p * 1000,
      trtf_norm  = t_normal["trtf"] * 1000,
      trtf_perm  = t_trtf_p * 1000,
      ks_norm    = t_normal["ks"] * 1000,
      ks_perm    = t_ks_p * 1000,
      ttm_norm   = t_normal["ttm"] * 1000,
      ttm_perm   = t_ttm_p * 1000
    )

  tab_all <- tab_all %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  header_lvl1 <- c(" " = 2,
                   "true" = 2,
                   "trtf" = 2,
                   "ks"   = 2,
                   "ttm"  = 2)
  header_lvl2 <- c(
    "dim" = 1,
    "distr" = 1,
    "normal" = 1, "permu" = 1,
    "normal" = 1, "permu" = 1,
    "normal" = 1, "permu" = 1,
    "normal" = 1, "permu" = 1
  )

  tab_all %>%
    kbl(caption = "Average logL", align = "c", booktabs = TRUE) %>%
    add_header_above(header_lvl1, align = "c") %>%
    add_header_above(header_lvl2, align = "c") %>%
    kable_styling(full_width = FALSE, position = "center") -> kbl_out

  attr(kbl_out, "tab_data") <- tab_all
  kbl_out
}


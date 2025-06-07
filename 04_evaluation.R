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

#' Summen-Zeile an Tabelle anhängen
#'
#' Fügt einer Log-Likelihood-Tabelle eine Zeile mit der
#' Gesamt-Summe der numerischen Spalten hinzu. Die Spalte
#' `dim` erhält den übergebenen Bezeichner.
#'
#' @param tab Datenrahmen mit numerischen Spalten
#' @param label Zeichenkette für die `dim`-Spalte
#' @return Datenrahmen mit zusätzlicher Summen-Zeile
#' @export
add_sum_row <- function(tab, label = "k") {
  stopifnot(is.data.frame(tab))
  sum_row <- setNames(vector("list", ncol(tab)), names(tab))
  for (nm in names(tab)) {
    if (nm == "dim") {
      sum_row[[nm]] <- label
    } else if (is.numeric(tab[[nm]])) {
      sum_row[[nm]] <- sum(tab[[nm]])
    } else {
      sum_row[[nm]] <- NA
    }
  }
  rbind(tab, as.data.frame(sum_row, stringsAsFactors = FALSE))
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
                                M_TRUE_p, M_TRTF_p, M_KS_p, M_TTM_p = NULL,
                                X_te_p, t_normal) {
  if (!requireNamespace("kableExtra", quietly = TRUE))
    install.packages("kableExtra", repos = "https://cloud.r-project.org")


  t_true_p <- system.time(logL_TRUE_dim(M_TRUE_p, X_te_p))["elapsed"]
  t_trtf_p <- system.time(logL_TRTF_dim(M_TRTF_p, X_te_p))["elapsed"]
  t_ks_p   <- system.time(logL_KS_dim(M_KS_p,   X_te_p))["elapsed"]
  if (!is.null(M_TTM_p)) {
    t_ttm_p <- system.time(logL_TTM_dim(M_TTM_p, X_te_p))["elapsed"]
  } else {
    t_ttm_p <- NA_real_
  }

  clean_cols <- function(df, suffix, has_ttm = TRUE) {
    df <- df %>%
      rename(
        distr = distribution,
        true  = logL_baseline,
        trtf  = logL_trtf,
        ks    = logL_ks
      )
    if (has_ttm && "logL_ttm" %in% names(df)) {
      df <- df %>% rename(ttm = logL_ttm)
      rn <- c("true", "trtf", "ks", "ttm")
    } else {
      rn <- c("true", "trtf", "ks")
    }
    df %>% rename_with(~ paste0(.x, "_", suffix), rn)
  }

  tab_norm <- clean_cols(tab_normal, "norm", !is.null(M_TTM_p))
  tab_perm <- clean_cols(tab_perm,   "perm", !is.null(M_TTM_p))

  tab_all <- left_join(tab_norm, tab_perm, by = c("dim", "distr"))
  tab_all <- tab_all %>%
    mutate(across(where(is.numeric), replace_na, 0)) %>%
    mutate(distr = replace_na(distr, "SUM"))

  tab_all <- tab_all %>%
    rename_with(~ sub("_norm$", "", .x), ends_with("_norm"))

  runtime_row <- list(
    dim       = "runtime (ms)",
    distr     = "",
    true      = t_normal["true"] * 1000,
    true_perm = t_true_p * 1000,
    trtf      = t_normal["trtf"] * 1000,
    trtf_perm = t_trtf_p * 1000,
    ks        = t_normal["ks"] * 1000,
    ks_perm   = t_ks_p * 1000
  )
  if (!is.null(M_TTM_p)) {
    runtime_row$ttm <- t_normal["ttm"] * 1000
    runtime_row$ttm_perm <- t_ttm_p * 1000
  }

  tab_all <- tab_all %>% add_row(!!!runtime_row)

  tab_all <- tab_all %>%
    identity()

  num_cols <- vapply(tab_all, is.numeric, logical(1))
  for (col in names(tab_all)[num_cols]) {
    tab_all[-nrow(tab_all), col] <- round(tab_all[-nrow(tab_all), col], 3)
    tab_all[nrow(tab_all), col] <- round(tab_all[nrow(tab_all), col], 0)
  }

  if (!is.null(M_TTM_p)) {
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
  } else {
    header_lvl1 <- c(" " = 2,
                     "true" = 2,
                     "trtf" = 2,
                     "ks"   = 2)
    header_lvl2 <- c(
      "dim" = 1,
      "distr" = 1,
      "normal" = 1, "permu" = 1,
      "normal" = 1, "permu" = 1,
      "normal" = 1, "permu" = 1
    )
  }

  tab_all %>%
    kbl(caption = "Average logL", align = "c", booktabs = TRUE) %>%
    add_header_above(header_lvl1, align = "c") %>%
    add_header_above(header_lvl2, align = "c") %>%
    kable_styling(full_width = FALSE, position = "center") -> kbl_out

  attr(kbl_out, "tab_data") <- tab_all
  kbl_out
}


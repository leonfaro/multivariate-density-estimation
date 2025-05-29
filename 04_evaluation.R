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

# Order learning heuristics for triangular transport maps (base R only)

#' Save permutation as CSV
#' @param perm integer permutation vector (length K)
#' @param path output CSV path (default 'artifacts/order_perm.csv')
save_perm <- function(perm, path = "artifacts/order_perm.csv") {
  stopifnot(is.integer(perm) || is.numeric(perm))
  perm <- as.integer(perm)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- data.frame(k = seq_along(perm), perm = perm)
  utils::write.csv(df, file = path, row.names = FALSE)
  invisible(path)
}

#' Learn a variable ordering via pivoted Cholesky
#'
#' @param S_or_X list with X_tr, X_te or numeric matrix X (train-only used)
#' @param seed integer RNG seed for determinism
#' @param method one of c("chol_pivot","identity")
#' @param gaussianize logical; if TRUE apply PIT->qnorm, else z-standardize only
#' @return list(perm=integer, method=character, gaussianize=logical)
learn_ordering <- function(S_or_X, seed = 42L, method = c("chol_pivot", "identity"), gaussianize = TRUE) {
  method <- match.arg(method)
  set.seed(as.integer(seed))

  # Extract training matrix
  if (is.list(S_or_X) && !is.null(S_or_X$X_tr)) {
    X_tr <- as.matrix(S_or_X$X_tr)
  } else {
    X_tr <- as.matrix(S_or_X)
  }
  stopifnot(is.matrix(X_tr), is.numeric(X_tr))
  N <- nrow(X_tr); K <- ncol(X_tr)
  if (K < 2L || method == "identity") {
    perm <- seq_len(K)
    save_perm(perm)
    return(list(perm = as.integer(perm), method = method, gaussianize = gaussianize))
  }

  # Train-only standardization
  mu <- colMeans(X_tr)
  sigma <- apply(X_tr, 2, sd) + .Machine$double.eps
  Xs <- sweep(sweep(X_tr, 2, mu, "-"), 2, sigma, "/")

  # Gaussianize via PIT with ties.average and clipping to (0,1)
  if (isTRUE(gaussianize)) {
    eps <- 1 / (N + 1)
    U <- matrix(NA_real_, N, K)
    for (k in seq_len(K)) {
      r <- rank(Xs[, k], ties.method = "average")
      u <- r / (N + 1)
      # Clip to [eps, 1-eps]
      u <- pmin(pmax(u, eps), 1 - eps)
      U[, k] <- u
    }
    Z <- stats::qnorm(U)
  } else {
    Z <- Xs
  }

  # Correlation and pivoted Cholesky
  R <- tryCatch(stats::cor(Z), error = function(e) NULL)
  perm <- seq_len(K)
  if (is.null(R) || any(!is.finite(R))) {
    warning("cor() failed or produced non-finite values; returning identity permutation")
  } else {
    cp <- tryCatch(chol(R, pivot = TRUE), error = function(e) NULL)
    if (is.null(cp)) {
      warning("chol(pivot=TRUE) failed; returning identity permutation")
    } else {
      pv <- attr(cp, "pivot")
      if (is.null(pv) || length(pv) != K) {
        warning("pivot attribute missing or invalid; returning identity permutation")
      } else {
        perm <- as.integer(pv)
      }
    }
  }

  save_perm(perm)
  list(perm = as.integer(perm), method = method, gaussianize = gaussianize)
}


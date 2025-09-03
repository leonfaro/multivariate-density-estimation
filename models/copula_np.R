# Nonparametric Copula Baseline (robust runner)
#
# Primary mode (2D, labeled): per-class univariate KDE marginals (dims 1 and 2)
# and a bivariate kernel copula on the empirical PIT (kde1d + kdecopula).
# Fallback mode (generic K-D, unlabeled or missing packages): independent
# kernel density product using base::density per dimension (no extra deps).

#' Fit per-class nonparametric copula (2D)
#'
#' @param S List with `X_tr` (numeric matrix with 2 columns) and `y_tr` (labels)
#' @param seed Integer RNG seed for determinism
#' @return An object of class 'copula_np' with fields `classes`, `by_class`, `seed`
#' @details
#' Per class y, the method:
#' - fits kde1d on each marginal X1,X2 of the class-specific subset (train only)
#' - computes pseudo-observations U via mid-ranks with ties.average, scaled by (n+1)
#' - clamps U to [1/(n+1), n/(n+1)]
#' - fits a kernel copula (kdecopula::kdecop) on U with method='TLL2'
#' The returned object stores for each class y: list(m1, m2, cop, n_tr, eps).
fit_copula_np <- function(S, seed = 42) {
  stopifnot(is.list(S))
  X_tr <- as.matrix(S$X_tr)
  y_tr <- S$y_tr
  K <- ncol(X_tr)
  has_pkgs <- (requireNamespace("kde1d", quietly = TRUE) && requireNamespace("kdecopula", quietly = TRUE))
  use_cop <- (K == 2L && !is.null(y_tr) && length(y_tr) == nrow(X_tr) && has_pkgs)

  set.seed(as.integer(seed))
  if (!use_cop) {
    # Fallback: product of univariate KDEs (base::density)
    dens <- vector("list", K)
    for (k in seq_len(K)) {
      dk <- stats::density(X_tr[, k], bw = "nrd0", n = 2048,
                           from = min(X_tr[, k]), to = max(X_tr[, k]))
      dens[[k]] <- list(x = dk$x, y = pmax(dk$y, .Machine$double.xmin))
    }
    return(structure(list(mode = "kde_product", dens = dens, seed = as.integer(seed)),
                    class = "copula_np"))
  }

  lbl <- sort(unique(as.vector(y_tr)))
  by_class <- setNames(vector("list", length(lbl)), as.character(lbl))

  # helper: mid-ranks scaled by (n+1)
  .midranks01 <- function(x) {
    n <- length(x)
    r <- rank(x, ties.method = "average")
    u <- r / (n + 1)
    eps <- 1 / (n + 1)
    pmin(pmax(u, eps), 1 - eps)
  }

  for (i in seq_along(lbl)) {
    y <- lbl[[i]]
    idx <- which(y_tr == y)
    Xi <- X_tr[idx, , drop = FALSE]
    n_i <- nrow(Xi)
    if (n_i < 5L) stop(sprintf("Class '%s' has too few samples (n=%d) to fit KDE/copula", as.character(y), n_i))
    eps <- 1 / (n_i + 1)

    # Fit univariate KDE marginals (store fit objects)
    m1 <- kde1d::kde1d(Xi[, 1])
    m2 <- kde1d::kde1d(Xi[, 2])

    # Empirical PIT via mid-ranks per class
    u1 <- .midranks01(Xi[, 1])
    u2 <- .midranks01(Xi[, 2])
    stopifnot(length(u1) == n_i, length(u2) == n_i,
              all(is.finite(u1)), all(is.finite(u2)))
    U <- cbind(u1, u2)

    # Fit kernel copula on U (TLL2)
    cop <- kdecopula::kdecop(U, method = "TLL2")

    # store per-class entry
    by_class[[as.character(y)]] <- list(
      m1 = m1,
      m2 = m2,
      cop = cop,
      n_tr = n_i,
      eps = eps
    )
  }

  # Basic numeric sanity for stored numeric fields
  for (k in names(by_class)) {
    ent <- by_class[[k]]
    if (!is.finite(ent$n_tr) || !is.finite(ent$eps)) {
      stop(sprintf("Non-finite numeric slot in by_class['%s']", k))
    }
  }

  structure(list(mode = "copula2d", classes = lbl, by_class = by_class, seed = as.integer(seed)),
            class = "copula_np")
}

# Optional: compact printer
print.copula_np <- function(x, ...) {
  cat(sprintf("copula_np model with %d classes; seed=%s\n",
              length(x$classes), as.character(x$seed)))
  cls <- as.character(x$classes)
  for (y in cls) {
    ent <- x$by_class[[as.character(y)]]
    cat(sprintf("  class %s: n_tr=%d, eps=%.4g\n", y, ent$n_tr, ent$eps))
  }
  invisible(x)
}

#' Predict log-density for copula_np model
#'
#' @param object copula_np model returned by fit_copula_np
#' @param newdata numeric matrix with 2 columns
#' @param y class labels (length must match nrow(newdata))
#' @param type one of 'logdensity_by_dim' or 'logdensity'
#' @return N x 2 matrix (by_dim) or length-N vector (joint)
predict.copula_np <- function(object, newdata, y, type = c("logdensity", "logdensity_by_dim")) {
  type <- match.arg(type)
  stopifnot(inherits(object, "copula_np"))
  X <- as.matrix(newdata)
  N <- nrow(X)
  if (!is.character(object$mode)) object$mode <- if (!is.null(object$by_class)) "copula2d" else "kde_product"

  if (identical(object$mode, "kde_product")) {
    # Ignore labels; per-dim KDE product
    K <- ncol(X)
    LD <- matrix(NA_real_, N, K)
    for (k in seq_len(K)) {
      d <- object$dens[[k]]
      fx <- stats::approx(d$x, d$y, xout = X[, k], rule = 2, ties = "ordered")$y
      LD[, k] <- log(pmax(fx, .Machine$double.xmin))
    }
    if (type == "logdensity_by_dim") return(LD)
    return(rowSums(LD))
  }

  # copula2d mode
  if (!(is.matrix(X) && is.numeric(X) && ncol(X) == 2L)) {
    stop("copula2d: newdata must be a numeric matrix with exactly 2 columns")
  }
  if (missing(y) || length(y) != N) {
    stop("copula2d: y must be provided and have the same length as nrow(newdata)")
  }
  if (!requireNamespace("kde1d", quietly = TRUE) || !requireNamespace("kdecopula", quietly = TRUE)) {
    stop("copula2d: packages 'kde1d' and 'kdecopula' are required")
  }
  LD <- matrix(0.0, N, 2L)
  yv <- as.vector(y)
  classes <- object$classes
  if (any(!(unique(yv) %in% classes))) stop("copula2d: label not in model$classes")
  tiny <- .Machine$double.xmin
  for (yy in classes) {
    id <- which(yv == yy)
    if (!length(id)) next
    comp <- object$by_class[[as.character(yy)]]
    Xi <- X[id, , drop = FALSE]
    f1 <- log(pmax(kde1d::dkde1d(Xi[, 1], comp$m1), tiny))
    F1 <- pmin(pmax(kde1d::pkde1d(Xi[, 1], comp$m1), comp$eps), 1 - comp$eps)
    f2 <- log(pmax(kde1d::dkde1d(Xi[, 2], comp$m2), tiny))
    F2 <- pmin(pmax(kde1d::pkde1d(Xi[, 2], comp$m2), comp$eps), 1 - comp$eps)
    lc <- log(pmax(kdecopula::dkdecop(cbind(F1, F2), comp$cop), tiny))
    LD[id, 1] <- f1 + 0.5 * lc
    LD[id, 2] <- f2 + 0.5 * lc
  }
  if (any(!is.finite(LD))) stop("NA/Inf in predict.copula_np output")
  if (type == "logdensity_by_dim") return(LD)
  rowSums(LD)
}

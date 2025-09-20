# Core helpers shared by TTM variants (marginal, separable)

# --------- Utilities and shared helpers ---------
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# Train-only standardization utilities
standardize_train_only <- function(Xtr, X = NULL) {
  stopifnot(is.matrix(Xtr) || is.data.frame(Xtr))
  Xtr <- as.matrix(Xtr)
  mu <- colMeans(Xtr)
  sigma <- apply(Xtr, 2, sd) + .Machine$double.eps
  if (is.null(X)) {
    Xs <- sweep(sweep(Xtr, 2, mu, "-"), 2, sigma, "/")
  } else {
    Xs <- sweep(sweep(as.matrix(X), 2, mu, "-"), 2, sigma, "/")
  }
  list(Xs = Xs, mu = mu, sigma = sigma)
}

.ttm_std_apply <- function(mu, sigma, X) {
  X1 <- sweep(X, 2, mu, "-")
  sweep(X1, 2, sigma, "/")
}

# Design/basis builders (formerly in ttm_bases.R)

# Replicate a row vector into an n-row matrix
rep_row <- function(vec, n) {
  v <- as.numeric(vec)
  matrix(rep(v, each = n), nrow = n, ncol = length(v))
}

# Gauss–Legendre nodes and weights on [0,1]
gauss_legendre_nodes <- function(Q) {
  stopifnot(is.numeric(Q), length(Q) == 1L, Q >= 1, Q == as.integer(Q))
  Q <- as.integer(Q)
  if (Q == 1L) return(list(nodes = 0.5, weights = 1.0))
  i <- seq_len(Q - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, Q, Q)
  for (k in i) {
    J[k, k + 1L] <- b[k]
    J[k + 1L, k] <- b[k]
  }
  e <- eigen(J, symmetric = TRUE)
  # Map from [-1,1] to [0,1]
  x <- (e$values + 1) / 2
  w <- (2 * (e$vectors[1L, ]^2)) / 2
  list(nodes = as.numeric(x), weights = as.numeric(w))
}
nodes_GL <- gauss_legendre_nodes

# Monotone 1-D basis in x: columns [x, erf(x)]
build_f <- function(x) {
  x <- as.numeric(x)
  cbind(x, 2 * stats::pnorm(x * sqrt(2)) - 1)
}

# Derivative of build_f: columns [1, 2/sqrt(pi)*exp(-x^2)]
d_build_f <- function(x) {
  x <- as.numeric(x)
  cbind(rep(1, length(x)), 2 / sqrt(pi) * exp(-x^2))
}

# Polynomial features in predecessors (without x_k)
build_g <- function(X_prev, deg = 3L, mask = NULL) {
  X_prev <- as.matrix(X_prev)
  stopifnot(deg == as.integer(deg), deg >= 0)
  N <- nrow(X_prev)
  p <- if (is.null(dim(X_prev))) 1L else ncol(X_prev)
  if (p == 0L) return(matrix(0, nrow = N, ncol = 0))
  out <- matrix(1, nrow = N, ncol = 1)
  if (deg >= 1L) {
    for (j in seq_len(p)) {
      xj <- X_prev[, j]
      for (d in seq_len(deg)) {
        out <- cbind(out, xj^d)
      }
    }
  }
  if (!is.null(mask)) {
    M <- try(as.matrix(mask), silent = TRUE)
    if (is.matrix(M)) {
      if (ncol(M) * nrow(M) == (ncol(out) - 1L)) {
        keep <- rep(TRUE, 1 + nrow(M) * ncol(M))
        idx <- 2L
        for (j in seq_len(nrow(M))) for (d in seq_len(ncol(M))) { if (isFALSE(M[j, d])) keep[idx] <- FALSE; idx <- idx + 1L }
        if (any(!keep)) out[, !keep, drop = FALSE] <- 0
      }
    }
  }
  out
}

# Row-wise Khatri–Rao (row-wise Kronecker)
.kr_rowwise <- function(A, B) {
  stopifnot(is.matrix(A), is.matrix(B), nrow(A) == nrow(B))
  N <- nrow(A); a <- ncol(A); b <- ncol(B)
  C <- matrix(0, nrow = N, ncol = a * b)
  for (j in seq_len(a)) {
    cols <- ((j - 1L) * b + 1L):(j * b)
    C[, cols] <- B * A[, j]
  }
  C
}

# Uniform open B-spline basis of degree 'degree' with 'df' basis functions
.bspline_basis_uniform <- function(t, df = 8L, degree = 3L, boundary = NULL) {
  stopifnot(df == as.integer(df), degree == as.integer(degree), df >= degree + 1L)
  t <- as.numeric(t)
  if (is.null(boundary)) boundary <- range(t)
  a <- boundary[1]; b <- boundary[2]
  n_basis <- as.integer(df); p <- as.integer(degree)
  n_interior <- n_basis - p - 1L
  if (n_interior <= 0L) interior <- numeric(0) else interior <- seq(a, b, length.out = n_interior + 2L)[-c(1L, n_interior + 2L)]
  knots <- c(rep(a, p + 1L), interior, rep(b, p + 1L))
  N0 <- function(x, i) as.numeric(knots[i] <= x & x < knots[i + 1L])
  N0b <- function(x, i) { y <- N0(x, i); if (i == (length(knots) - 1L)) y[x >= knots[length(knots)]] <- 1; y }
  M <- length(t); n <- n_basis
  B <- matrix(0, nrow = M, ncol = n)
  for (i in seq_len(n)) B[, i] <- N0b(t, i)
  if (p >= 1L) {
    for (q in 1:p) {
      Bq <- matrix(0, nrow = M, ncol = n)
      for (i in seq_len(n)) {
        denom1 <- knots[i + q] - knots[i]
        denom2 <- knots[i + q + 1L] - knots[i + 1L]
        a1 <- if (denom1 > 0) (t - knots[i]) / denom1 else 0
        a2 <- if (denom2 > 0) (knots[i + q + 1L] - t) / denom2 else 0
        Bq[, i] <- a1 * B[, i] + a2 * if (i + 1L <= n) B[, i + 1L] else 0
      }
      B <- Bq
    }
  }
  B
}

# --------- Ordering heuristics (formerly order_heuristics.R) ---------
save_perm <- function(perm, path = "artifacts/order_perm.csv") {
  stopifnot(is.integer(perm) || is.numeric(perm))
  perm <- as.integer(perm)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- data.frame(k = seq_along(perm), perm = perm)
  utils::write.csv(df, file = path, row.names = FALSE)
  invisible(path)
}

learn_ordering <- function(S_or_X, seed = 42L, method = c("chol_pivot", "identity"), gaussianize = TRUE) {
  method <- match.arg(method)
  set.seed(as.integer(seed))
  if (is.list(S_or_X) && !is.null(S_or_X$X_tr)) X_tr <- as.matrix(S_or_X$X_tr) else X_tr <- as.matrix(S_or_X)
  stopifnot(is.matrix(X_tr), is.numeric(X_tr))
  N <- nrow(X_tr); K <- ncol(X_tr)
  if (K < 2L || method == "identity") {
    perm <- seq_len(K); save_perm(perm); return(list(perm = as.integer(perm), method = method, gaussianize = gaussianize))
  }
  mu <- colMeans(X_tr); sigma <- apply(X_tr, 2, sd) + .Machine$double.eps
  Xs <- sweep(sweep(X_tr, 2, mu, "-"), 2, sigma, "/")
  if (isTRUE(gaussianize)) {
    eps <- 1 / (N + 1)
    U <- matrix(NA_real_, N, K)
    for (k in seq_len(K)) {
      r <- rank(Xs[, k], ties.method = "average"); u <- r / (N + 1)
      U[, k] <- pmin(pmax(u, eps), 1 - eps)
    }
    Z <- stats::qnorm(U)
  } else Z <- Xs
  R <- tryCatch(stats::cor(Z), error = function(e) NULL)
  perm <- seq_len(K)
  if (!is.null(R) && all(is.finite(R))) {
    cp <- tryCatch(chol(R, pivot = TRUE), error = function(e) NULL)
    if (!is.null(cp)) {
      pv <- attr(cp, "pivot"); if (!is.null(pv) && length(pv) == K) perm <- as.integer(pv)
    }
  }
  save_perm(perm)
  list(perm = as.integer(perm), method = method, gaussianize = gaussianize)
}

# --------- Algo dispatch and forward/Jacobian ---------

.ttm_algo_of <- function(model) {
  if (!is.null(model$algo)) return(model$algo)
  cl <- class(model)
  if ("ttm_marginal2" %in% cl || "ttm_marginal" %in% cl) return("marginal")
  if ("ttm_separable" %in% cl) return("separable")
  stop("Unknown TTM model class for core forward computation")
}

# ttm_forward: returns standardized forward map outputs and standardized Jacobian diagonal
ttm_forward <- function(model, X) {
  stopifnot(is.matrix(X))
  mu <- model$mu; sigma <- model$sigma
  stopifnot(is.numeric(mu), is.numeric(sigma), length(mu) == ncol(X), length(sigma) == ncol(X))
  Xs <- .ttm_std_apply(mu, sigma, X)
  N <- nrow(Xs); K <- ncol(Xs)
  algo <- .ttm_algo_of(model)
  Z <- matrix(0, N, K)
  J <- matrix(0, N, K)
  if (identical(algo, "marginal")) {
    if (!is.null(model$coeffs)) {
      for (k in seq_len(K)) {
        a <- as.numeric(model$coeffs[[k]]["a"])
        b <- as.numeric(model$coeffs[[k]]["b"])
        Z[, k] <- a + b * Xs[, k]
        J[, k] <- b
      }
    } else if (!is.null(model$coeffA) && !is.null(model$coeffB)) {
      b <- exp(model$coeffA); a <- model$coeffB
      for (k in seq_len(K)) { Z[, k] <- a[k] + b[k] * Xs[, k]; J[, k] <- b[k] }
    } else stop("Marginal model missing required coefficients")
  } else if (identical(algo, "separable")) {
    deg <- if (!is.null(model$degree_g)) model$degree_g else 2L
    for (k in seq_len(K)) {
      x_prev <- if (k > 1) Xs[, 1:(k-1), drop = FALSE] else matrix(0, N, 0)
      xk <- Xs[, k]
      if (k > 1) {
        P_non <- build_g(x_prev, deg)
      } else {
        if (!is.null(model$degree_g) && model$degree_g == 0L && !is.null(model$coeffs[[k]]$c_non)) {
          P_non <- matrix(1, nrow = N, ncol = 1L)
        } else {
          P_non <- matrix(0, N, 0)
        }
      }
      P_mon <- build_f(xk)
      B <- d_build_f(xk)
      if (!is.null(model$coeffs[[k]]$c_mon) && length(model$coeffs[[k]]$c_mon) < ncol(P_mon)) {
        m_keep <- length(model$coeffs[[k]]$c_mon)
        P_mon <- P_mon[, seq_len(m_keep), drop = FALSE]
        B <- B[, seq_len(m_keep), drop = FALSE]
      }
      c_non <- model$coeffs[[k]]$c_non
      c_mon <- model$coeffs[[k]]$c_mon
      gk <- if (ncol(P_non) > 0) as.numeric(P_non %*% c_non) else rep(0, N)
      fk <- as.numeric(P_mon %*% c_mon)
      Z[, k] <- gk + fk
      J[, k] <- as.numeric(B %*% c_mon)
    }
  } else stop("Unsupported algo")
  list(Z = Z, J = J)
}

# Per-dimension log-density under N(0,1)^K reference and train-only standardization
ttm_ld_by_dim <- function(model, X) {
  out <- ttm_forward(model, X)
  Z <- out$Z; J <- out$J
  stopifnot(is.matrix(Z), is.matrix(J), all(dim(Z) == dim(J)))
  N <- nrow(Z); K <- ncol(Z)
  C <- -0.5 * log(2 * pi)
  LJ <- log(J) - matrix(log(model$sigma), nrow = N, ncol = K, byrow = TRUE)
  (-0.5) * (Z^2) + C + LJ
}

# Unified design/basis builders for TTM variants (pure base R)

# Robust repo-local source helper (works from tests or scripts)
.source_repo_file <- function(rel_path) {
  cand <- character(0)
  # 1) relative to current working directory
  cand <- c(cand, rel_path)
  # 2) relative to this file if available
  this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
  if (is.character(this_file) && nzchar(this_file) && !is.na(this_file)) {
    cand <- c(cand, file.path(dirname(this_file), basename(rel_path)))
  }
  # 3) via root_path if defined
  if (exists("root_path")) cand <- c(cand, file.path(root_path, rel_path))
  # 4) walk up to find 00_globals.R as repo root
  p <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:8) {
    if (file.exists(file.path(p, "00_globals.R"))) { cand <- c(cand, file.path(p, rel_path)); break }
    np <- dirname(p); if (identical(np, p)) break; p <- np
  }
  cand <- unique(cand)
  for (pth in cand) { if (file.exists(pth)) { source(pth); return(invisible(TRUE)) } }
  invisible(FALSE)
}

# Helper: replicate a row vector into an n-row matrix
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

# Alias commonly used name
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
# Returns intercept + powers up to degree 'deg' for each column in X_prev
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
  # Optional masking of polynomial powers per predecessor (p x deg logical)
  if (!is.null(mask)) {
    M <- try(as.matrix(mask), silent = TRUE)
    if (is.matrix(M) && nrow(M) == p && ncol(M) == deg) {
      # Columns layout after intercept: for j=1..p, blocks of size 'deg' in order d=1..deg
      keep <- rep(TRUE, 1 + p * deg)
      idx <- 2L
      for (j in seq_len(p)) {
        for (d in seq_len(deg)) {
          if (isFALSE(M[j, d])) keep[idx] <- FALSE
          idx <- idx + 1L
        }
      }
      # Zero-out excluded columns to keep API/stats consistent (avoid changing dims)
      if (any(!keep)) {
        out[, !keep, drop = FALSE] <- 0
      }
    }
  }
  out
}

# Row-wise Khatri–Rao (row-wise Kronecker): for each row i, vec(A[i,] ⊗ B[i,])
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
  # Open uniform knot vector
  n_basis <- as.integer(df)
  p <- as.integer(degree)
  n_knots <- n_basis + p + 1L
  n_interior <- n_basis - p - 1L
  if (n_interior <= 0L) {
    interior <- numeric(0)
  } else {
    interior <- seq(a, b, length.out = n_interior + 2L)[-c(1L, n_interior + 2L)]
  }
  knots <- c(rep(a, p + 1L), interior, rep(b, p + 1L))
  # Cox–de Boor recursion
  N0 <- function(x, i) as.numeric(knots[i] <= x & x < knots[i + 1L])
  # Handle x==b to be included in last basis
  N0b <- function(x, i) {
    y <- N0(x, i)
    if (i == (length(knots) - 1L)) y[x >= knots[length(knots)]] <- 1
    y
  }
  # Initialize basis for degree 0
  M <- length(t)
  n <- n_basis
  B <- matrix(0, nrow = M, ncol = n)
  for (i in seq_len(n)) {
    B[, i] <- N0b(t, i)
  }
  # Higher degrees
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

if (!exists(".bspline_basis_open_knots") || !exists("build_symmetric_knots")) {
  .source_repo_file(file.path("models", "ttm", "crossterm_basis.R"))
}
if (!exists("build_t_rbf")) {
  .source_repo_file(file.path("models", "ttm", "rbf_basis.R"))
}

# build_h: tensor basis between B-splines in t and polynomial features in X_prev
build_h <- function(t, X_prev, spec = list(df = 8L, degree = 3L, deg_g = 3L)) {
  t <- as.numeric(t)
  X_prev <- as.matrix(X_prev)
  N <- length(t)
  stopifnot(nrow(X_prev) == N || ncol(X_prev) == 0L)
  df <- if (!is.null(spec$df)) as.integer(spec$df) else 8L
  degree <- if (!is.null(spec$degree)) as.integer(spec$degree) else 3L
  deg_g <- if (!is.null(spec$deg_g)) as.integer(spec$deg_g) else 3L
  # Basis in t: switch on spec$kind; default is B-spline
  kind <- if (!is.null(spec$kind)) as.character(spec$kind) else "bspline"
  if (identical(kind, "rbf") && exists("build_t_rbf")) {
    # RBF integrated basis with smooth tail-linear behavior internally
    Bt <- build_t_rbf(t, df = df, prob_hi = 0.999)
  } else {
    # Support symmetric knots and constant tails if provided in spec
    if (!is.null(spec$knots_t) || !is.null(spec$boundary_t)) {
      interior <- spec$knots_t %||% numeric(0)
      boundary <- spec$boundary_t %||% range(t)
      Bt <- .bspline_basis_open_knots(t, interior = interior, degree = degree, boundary = boundary)
      if (isTRUE(spec$tail_const %||% TRUE)) {
        Tfeat <- tail_plateau_features(t, boundary)
        Bt <- cbind(Bt, Tfeat)
      }
    } else {
      Bt <- .bspline_basis_uniform(t, df = df, degree = degree)
    }
  }
  Gx <- if (ncol(X_prev) > 0L) build_g(X_prev, deg = deg_g) else matrix(1, nrow = N, ncol = 1L)
  .kr_rowwise(Bt, Gx)
}

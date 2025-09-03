## Cross-term t-basis helpers with symmetric knots and tail features

# Open B-spline basis from explicit interior knots and boundary [a,b]
.bspline_basis_open_knots <- function(t, interior = numeric(0), degree = 3L, boundary = NULL) {
  stopifnot(degree == as.integer(degree), degree >= 0)
  t <- as.numeric(t)
  if (is.null(boundary)) {
    if (length(interior) == 0L) {
      boundary <- range(t)
    } else {
      boundary <- range(c(interior, t))
    }
  }
  a <- boundary[1]; b <- boundary[2]
  # Ensure strict ordering and within bounds
  interior <- sort(as.numeric(interior))
  interior <- interior[interior > a & interior < b]
  p <- as.integer(degree)
  n_basis <- length(interior) + p + 1L
  # Open uniform knot vector induced by interior+boundary
  knots <- c(rep(a, p + 1L), interior, rep(b, p + 1L))
  # degree-0 basis
  M <- length(t)
  N0 <- function(x, i) as.numeric(knots[i] <= x & x < knots[i + 1L])
  N0b <- function(x, i) {
    y <- N0(x, i)
    if (i == (length(knots) - 1L)) y[x >= knots[length(knots)]] <- 1
    y
  }
  B <- matrix(0, nrow = M, ncol = n_basis)
  for (i in seq_len(n_basis)) B[, i] <- N0b(t, i)
  if (p >= 1L) {
    for (q in 1:p) {
      Bq <- matrix(0, nrow = M, ncol = n_basis)
      for (i in seq_len(n_basis)) {
        denom1 <- knots[i + q] - knots[i]
        denom2 <- knots[i + q + 1L] - knots[i + 1L]
        a1 <- if (denom1 > 0) (t - knots[i]) / denom1 else 0
        a2 <- if (denom2 > 0) (knots[i + q + 1L] - t) / denom2 else 0
        Bq[, i] <- a1 * B[, i] + a2 * if (i + 1L <= n_basis) B[, i + 1L] else 0
      }
      B <- Bq
    }
  }
  attr(B, "knots") <- knots
  attr(B, "boundary") <- c(a, b)
  B
}

# Construct symmetric interior knots about 0 from empirical t-samples
build_symmetric_knots <- function(t_sample, df = 8L, degree = 3L, prob_hi = 0.995) {
  t_sample <- as.numeric(t_sample)
  t_sample <- t_sample[is.finite(t_sample)]
  if (length(t_sample) == 0L) return(list(interior = numeric(0), boundary = c(-1, 1)))
  p <- as.integer(degree)
  n_basis <- as.integer(df)
  n_int <- max(0L, n_basis - p - 1L)
  T <- as.numeric(stats::quantile(abs(t_sample), probs = prob_hi, names = FALSE))
  T <- max(T, 1.0)
  if (n_int == 0L) return(list(interior = numeric(0), boundary = c(-T, T)))
  m_half <- n_int %/% 2L
  probs <- if (m_half > 0L) seq(1, m_half) / (m_half + 1) else numeric(0)
  pos <- if (length(probs) > 0) as.numeric(stats::quantile(abs(t_sample), probs = probs, names = FALSE)) else numeric(0)
  pos <- pmin(pmax(pos, 1e-8), T * 0.999)
  if ((n_int %% 2L) == 1L) {
    interior <- c(-rev(pos), 0, pos)
  } else {
    interior <- c(-rev(pos), pos)
  }
  list(interior = as.numeric(interior), boundary = c(-T, T))
}

# Tail features that produce constant h in tails (per side)
tail_plateau_features <- function(t, boundary) {
  a <- boundary[1]; b <- boundary[2]
  left  <- as.numeric(t <= a)
  right <- as.numeric(t >= b)
  cbind(left, right)
}


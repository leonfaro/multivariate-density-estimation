# Radial Basis helpers for t-basis (base R only)

# Local coordinates (dimensionless)
rbf_local_coords <- function(x, mu, sigma) {
  s <- pmax(as.numeric(sigma), .Machine$double.eps)
  (as.numeric(x) - as.numeric(mu)) / s
}

# Gaussian RBF in local coords
rbf <- function(xloc) {
  exp(-0.5 * (as.numeric(xloc)^2))
}

# Integral of Gaussian RBF in local coords, such that d/dxloc irbf = rbf
# Using relation erf(x/sqrt(2)) = 2*pnorm(x) - 1
irbf <- function(xloc) {
  sqrt(pi / 2) * (2 * stats::pnorm(as.numeric(xloc)) - 1)
}

# Left-end tail term (smooth, asymptotically linear with slope 1 as x -> -inf)
# xloc = (x - a) / sigma  (a is left boundary), tau = sigma
let_term <- function(xloc, sigma) {
  tau <- pmax(as.numeric(sigma), .Machine$double.eps)
  # LET(x) = -tau * log1p(exp(-xloc)) + tau * log(2)
  -tau * log1p(exp(-as.numeric(xloc))) + tau * log(2)
}

# Right-end tail term (smooth, asymptotically linear with slope 1 as x -> +inf)
# xloc = (x - b) / sigma  (b is right boundary), tau = sigma
ret_term <- function(xloc, sigma) {
  tau <- pmax(as.numeric(sigma), .Machine$double.eps)
  # RET(x) =  tau * log1p(exp( xloc)) - tau * log(2)
  tau * log1p(exp(as.numeric(xloc))) - tau * log(2)
}

# Place RBF knots by empirical quantiles; local scales by neighbor distances
# Returns list(mu = centers, sigma = local_scales)
place_rbf_knots <- function(xk, m) {
  xk <- as.numeric(xk)
  m <- as.integer(m)
  stopifnot(m >= 1L)
  if (length(xk) == 0L) return(list(mu = numeric(0), sigma = numeric(0)))
  # Probabilities avoid exact 0/1 to keep centers within data range
  probs <- if (m == 1L) 0.5 else seq(1, m) / (m + 1)
  mu <- as.numeric(stats::quantile(xk, probs = probs, names = FALSE, type = 7))
  mu <- sort(unique(mu))
  L <- length(mu)
  if (L == 0L) return(list(mu = numeric(0), sigma = numeric(0)))
  # Local sigma from neighbor distances (half of adjacent gaps)
  sig <- rep(NA_real_, L)
  if (L == 1L) {
    r <- stats::IQR(xk)
    if (!is.finite(r) || r <= 0) r <- sd(xk)
    if (!is.finite(r) || r <= 0) r <- max(1.0, diff(range(xk)))
    sig[1] <- r / 6
  } else {
    gaps <- diff(mu)
    sig[1] <- gaps[1] / 2
    sig[L] <- gaps[L - 1] / 2
    if (L > 2L) for (j in 2:(L - 1)) sig[j] <- 0.5 * (mu[j + 1] - mu[j - 1]) / 2
  }
  sig <- pmax(sig, .Machine$double.eps)
  list(mu = mu, sigma = sig)
}

# Build t-basis with integrated RBF columns + smooth tail-linear terms
# df: number of RBF centers; prob_hi: sets tail boundaries via quantiles
build_t_rbf <- function(t, df = 8L, prob_hi = 0.995) {
  t <- as.numeric(t)
  N <- length(t)
  m <- as.integer(df)
  # RBF centers/scales
  knots <- place_rbf_knots(t, m)
  mu <- knots$mu; sg <- knots$sigma
  M <- length(mu)
  # iRBF columns (scale back to t via sigma_j)
  Phi <- if (M > 0L) {
    cols <- lapply(seq_len(M), function(j) {
      xloc <- rbf_local_coords(t, mu[j], sg[j])
      sg[j] * irbf(xloc)
    })
    do.call(cbind, cols)
  } else matrix(0, N, 0)
  # Tail boundaries a,b from quantiles
  a <- as.numeric(stats::quantile(t, probs = 1 - prob_hi, names = FALSE))
  b <- as.numeric(stats::quantile(t, probs = prob_hi, names = FALSE))
  if (!is.finite(a)) a <- min(t)
  if (!is.finite(b)) b <- max(t)
  if (a > b) { tmp <- a; a <- b; b <- tmp }
  # Tail smoothing scale from global spread
  tau <- max(.Machine$double.eps, 0.5 * stats::IQR(t))
  if (!is.finite(tau) || tau <= 0) tau <- max(.Machine$double.eps, sd(t))
  if (!is.finite(tau) || tau <= 0) tau <- 1.0
  # LET/RET terms
  xloc_L <- (t - a) / tau
  xloc_R <- (t - b) / tau
  LET <- let_term(xloc_L, tau)
  RET <- ret_term(xloc_R, tau)
  # Assemble basis: [iRBF..., LET, RET]
  out <- cbind(Phi, LET, RET)
  colnames(out) <- c(if (M > 0L) paste0("iRBF_", seq_len(M)) else character(0), "LET", "RET")
  stopifnot(is.matrix(out), nrow(out) == N, all(is.finite(out)))
  out
}


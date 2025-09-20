# True log-density for the two-moons generator (mixture of two arcs convolved with N(0, sigma^2 I))

gauss_legendre_0_pi <- function(Q) {
  if (Q <= 0L || Q != as.integer(Q)) stop("Q must be positive integer")
  if (Q == 1L) return(list(t = pi/2, w = pi))
  i <- seq_len(Q - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, Q, Q)
  for (k in i) { J[k, k + 1] <- b[k]; J[k + 1, k] <- b[k] }
  e <- eigen(J, symmetric = TRUE)
  u <- (e$values + 1) / 2  # in [0,1]
  t <- pi * u
  w <- (2 * (e$vectors[1, ]^2)) / 2  # on [0,1]
  w <- pi * w  # scale to [0,pi]
  list(t = as.numeric(t), w = as.numeric(w))
}

true_logdensity <- function(X, S, Q = 32L) {
  stopifnot(is.matrix(X), ncol(X) == 2)
  sigma <- if (!is.null(S$meta$sigma)) S$meta$sigma else if (!is.null(S$meta$noise)) S$meta$noise else 0.15
  gap <- if (!is.null(S$meta$gap)) S$meta$gap else 0.5
  R <- if (!is.null(S$meta$radius)) S$meta$radius else 1.0
  gl <- gauss_legendre_0_pi(Q)
  t <- gl$t; w <- gl$w
  # Curve means
  m1 <- cbind(R * cos(t), R * sin(t))
  m2 <- cbind(R * (1 - cos(t)), -R * sin(t) + gap)
  inv2s2 <- 1 / (2 * sigma^2)
  log_norm_const <- -log(2 * pi * sigma^2)
  N <- nrow(X)
  # Compute log φσ for each point against all nodes for both arcs
  # Efficient by broadcasting: for each q, compute distances to all X
  logp1 <- numeric(N)
  logp2 <- numeric(N)
  for (q in seq_along(t)) {
    d1 <- rowSums((X - matrix(m1[q, ], nrow = N, ncol = 2, byrow = TRUE))^2)
    d2 <- rowSums((X - matrix(m2[q, ], nrow = N, ncol = 2, byrow = TRUE))^2)
    lp1_q <- log(w[q]/pi) + log_norm_const - inv2s2 * d1
    lp2_q <- log(w[q]/pi) + log_norm_const - inv2s2 * d2
    if (q == 1L) {
      logp1 <- lp1_q; logp2 <- lp2_q
    } else {
      m1q <- pmax(logp1, lp1_q); logp1 <- m1q + log(exp(logp1 - m1q) + exp(lp1_q - m1q))
      m2q <- pmax(logp2, lp2_q); logp2 <- m2q + log(exp(logp2 - m2q) + exp(lp2_q - m2q))
    }
  }
  # Mix the two arcs with weight 0.5 each
  a <- log(0.5) + logp1
  b <- log(0.5) + logp2
  m <- pmax(a, b)
  joint <- m + log(exp(a - m) + exp(b - m))
  by_dim <- cbind(joint/2, joint/2)
  list(joint = joint, by_dim = by_dim)
}

# Conditional oracle log-density p(x | y) for two-moons
# Returns list(joint, by_dim) with by_dim columns summing to joint (within FP tolerance)
true_logdensity_conditional <- function(X, y, S, Q = 32L) {
  stopifnot(is.matrix(X), ncol(X) == 2, length(y) == nrow(X))
  sigma <- if (!is.null(S$meta$sigma)) S$meta$sigma else if (!is.null(S$meta$noise)) S$meta$noise else 0.15
  gap <- if (!is.null(S$meta$gap)) S$meta$gap else 0.5
  R <- if (!is.null(S$meta$radius)) S$meta$radius else 1.0
  gl <- gauss_legendre_0_pi(Q)
  t <- gl$t; w <- gl$w
  # Curve means for parameter t
  m1 <- cbind(R * cos(t), R * sin(t))
  m2 <- cbind(R * (1 - cos(t)), -R * sin(t) + gap)
  inv2s2 <- 1 / (2 * sigma^2)
  log_norm_const <- -log(2 * pi * sigma^2)
  N <- nrow(X)
  y <- as.integer(y)
  if (!all(y %in% c(1L, 2L))) stop("y must be integer labels in {1,2}")
  # Accumulate log-sum-exp across quadrature nodes for the selected arc per sample
  joint <- rep(-Inf, N)
  for (q in seq_along(t)) {
    # distances to means for both arcs at node q
    d1 <- rowSums((X - matrix(m1[q, ], nrow = N, ncol = 2, byrow = TRUE))^2)
    d2 <- rowSums((X - matrix(m2[q, ], nrow = N, ncol = 2, byrow = TRUE))^2)
    lp1_q <- log(w[q]) - log(pi) + log_norm_const - inv2s2 * d1
    lp2_q <- log(w[q]) - log(pi) + log_norm_const - inv2s2 * d2
    lp_q <- ifelse(y == 1L, lp1_q, lp2_q)
    m <- pmax(joint, lp_q)
    joint <- m + log(exp(joint - m) + exp(lp_q - m))
  }
  by_dim <- cbind(joint/2, joint/2)
  list(joint = joint, by_dim = by_dim)
}

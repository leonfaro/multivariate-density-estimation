# Per-k training driver for TTM variants (maps-from-samples, Eq. 38/39)
# Returns parameter list per dimension and meta.

source(file.path("R", "ttm_bases.R"))
source(file.path("R", "ttm_core.R"))

.opt_marginal_k <- function(xk, lambda = 0.0, eps = 1e-6) {
  # Closed-form matching marginal fit implementation (uses sample variance)
  v <- stats::var(xk) + 1e-12
  b <- max(eps, 1 / sqrt(v))
  a <- -b * mean(xk)
  list(a = a, b = b)
}

.opt_separable_k <- function(x_prev, xk, degree_g = 2L, lambda = 0.0, eps = 1e-6) {
  N <- length(xk)
  Phi_non <- if (ncol(x_prev) > 0) build_g(x_prev, deg = degree_g) else if (degree_g == 0L) matrix(1, N, 1L) else matrix(0, N, 0)
  if (degree_g == 0L) {
    v <- stats::var(xk) + 1e-12
    c_mon <- max(eps, 1 / sqrt(v))
    c_non <- matrix(-mean(xk) * c_mon, ncol = 1L)
    return(list(c_non = c_non, c_mon = c_mon))
  }
  Phi_mon <- build_f(xk)
  B <- d_build_f(xk)
  m_non <- ncol(Phi_non); m_mon <- ncol(Phi_mon); N <- nrow(Phi_mon)
  if (m_non > 0) {
    M <- solve(crossprod(Phi_non) + lambda * diag(m_non), t(Phi_non))
    A <- (diag(N) - Phi_non %*% M) %*% Phi_mon
    D <- M %*% Phi_mon
  } else {
    M <- matrix(0, 0, N)
    A <- Phi_mon
    D <- matrix(0, 0, m_mon)
  }
  fn <- function(c) {
    r <- A %*% c
    Bc <- B %*% c
    if (any(Bc <= 0)) return(Inf)
    q <- D %*% c
    0.5 * sum(r^2) - sum(log(Bc)) + 0.5 * lambda * (sum(q^2) + sum(c^2))
  }
  gr <- function(c) {
    r <- A %*% c
    Bc <- B %*% c
    q <- D %*% c
    as.numeric(t(A) %*% r - t(B) %*% (1 / Bc) + lambda * (t(D) %*% q + c))
  }
  c0 <- rep(1, m_mon)
  opt <- optim(c0, fn, gr, method = "L-BFGS-B", lower = rep(eps, m_mon))
  c_mon <- opt$par
  c_non <- if (m_non > 0) -M %*% (Phi_mon %*% c_mon) else numeric(0)
  list(c_non = c_non, c_mon = c_mon)
}

.opt_crossterm_k <- function(x_prev, xk, deg_g = 2L, df_t = 6L, lambda = 1e-3, Q = 16L, Hmax = 20) {
  N <- length(xk)
  Gx <- if (ncol(x_prev) > 0) build_g(x_prev, deg = deg_g) else matrix(1, N, 1L)
  m_g <- ncol(Gx)
  spec_h <- list(df = as.integer(df_t), degree = 3L, deg_g = as.integer(deg_g))
  gl <- gauss_legendre_nodes(Q)
  nodes <- gl$nodes; weights <- gl$weights
  Hstar <- build_h(xk, x_prev, spec_h)
  m_h <- ncol(Hstar)
  fn <- function(theta) {
    a <- if (m_g > 0) theta[seq_len(m_g)] else numeric(0)
    b <- theta[(m_g + 1L):length(theta)]
    g <- if (m_g > 0) as.numeric(Gx %*% a) else rep(0, N)
    I <- rep(0, N)
    for (q in seq_along(nodes)) {
      tq <- xk * nodes[q]
      Hq <- build_h(tq, x_prev, spec_h)
      v <- as.numeric(Hq %*% b)
      v <- pmax(pmin(v, Hmax), -Hmax)
      I <- I + weights[q] * exp(v)
    }
    I <- sign(xk) * abs(xk) * I
    S <- g + I
    h_star <- as.numeric(Hstar %*% b)
    h_star <- pmax(pmin(h_star, Hmax), -Hmax)
    0.5 * sum(S^2) - sum(h_star) + 0.5 * lambda * sum(theta^2)
  }
  theta0 <- c(rep(0, m_g), rep(0, m_h))
  opt <- optim(theta0, fn, method = "L-BFGS-B")
  a_hat <- if (m_g > 0) opt$par[seq_len(m_g)] else numeric(0)
  b_hat <- opt$par[(m_g + 1L):length(opt$par)]
  list(alpha = a_hat, beta = b_hat, spec_h = spec_h, nodes = nodes, weights = weights, Hmax = Hmax)
}

# Main training driver (per-k, pure loop)
train_per_k <- function(algo = c("marginal","separable","crossterm"), X_std, lambda = 0.0, deg = 2L, Q = 16L, seed = 42, Hmax = 20, mu = NULL, sigma = NULL) {
  set.seed(seed)
  algo <- match.arg(algo)
  stopifnot(is.matrix(X_std))
  N <- nrow(X_std); K <- ncol(X_std)
  params <- vector("list", K)
  t_train <- system.time({
    for (k in seq_len(K)) {
      x_prev <- if (k > 1) X_std[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
      xk <- X_std[, k]
      if (algo == "marginal") params[[k]] <- .opt_marginal_k(xk, lambda)
      if (algo == "separable") params[[k]] <- .opt_separable_k(x_prev, xk, degree_g = deg, lambda = lambda)
      if (algo == "crossterm") params[[k]] <- .opt_crossterm_k(x_prev, xk, deg_g = deg, df_t = 6L, lambda = lambda, Q = Q, Hmax = Hmax)
    }
  })[["elapsed"]]
  meta <- list(algo = algo,
               mu = if (is.null(mu)) rep(0, K) else mu,
               sigma = if (is.null(sigma)) rep(1, K) else sigma,
               time_train = t_train,
               lambda = lambda, deg = deg, Q = Q, seed = seed)
  list(params = params, meta = meta)
}

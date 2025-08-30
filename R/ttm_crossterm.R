# Cross-term triangular map: S_k(x) = g_k(x_<k>) + ∫_0^{x_k} exp(h_k(t, x_<k})) dt
# Minimal maps-from-samples fit using unified bases and core predictor.

source(file.path("R", "ttm_bases.R"))
source(file.path("R", "ttm_core.R"))

.std_stats <- function(X) {
  mu <- colMeans(X)
  sigma <- apply(X, 2, sd) + .Machine$double.eps
  list(mu = mu, sigma = sigma)
}

.standardize <- function(X, mu, sigma) sweep(sweep(X, 2, mu, "-"), 2, sigma, "/")

fit_ttm_crossterm <- function(data, deg_g = 2L, df_t = 6L, lambda = 1e-3, Q = 16L, Hmax = 20, maxit = 50L, seed = 42) {
  set.seed(seed)
  S_in <- if (is.list(data) && !is.null(data$X_tr) && !is.null(data$X_te)) data else {
    stopifnot(is.matrix(data))
    if (!exists("split_data")) source("02_split.R")
    list(X_tr = split_data(data, seed)$X_tr, X_te = split_data(data, seed)$X_te)
  }
  X_tr <- S_in$X_tr; X_te <- S_in$X_te
  N <- nrow(X_tr); K <- ncol(X_tr)

  time_train <- system.time({
    st <- .std_stats(X_tr); mu <- st$mu; sigma <- st$sigma
    Xs <- .standardize(X_tr, mu, sigma)
    gl <- gauss_legendre_nodes(Q)
    nodes <- gl$nodes; weights <- gl$weights
    spec_h <- list(df = as.integer(df_t), degree = 3L, deg_g = as.integer(deg_g))
    coeffs <- vector("list", K)
    for (k in seq_len(K)) {
      Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
      xk <- Xs[, k]
      Gx <- if (ncol(Xprev) > 0) build_g(Xprev, deg = deg_g) else matrix(1, N, 1L)
      m_g <- ncol(Gx)
      # h-basis uses the same poly degree in predecessors
      Hstar <- build_h(xk, Xprev, spec_h)  # N x (df * m_g_h)
      m_h <- ncol(Hstar)
      # Objective using numeric gradient (L-BFGS-B on fn only)
      fn_k <- function(theta) {
        a <- if (m_g > 0) theta[seq_len(m_g)] else numeric(0)
        b <- theta[(m_g + 1L):length(theta)]
        g <- if (m_g > 0) as.numeric(Gx %*% a) else rep(0, N)
        # Integral via GL on [0, xk]
        I <- rep(0, N)
        for (q in seq_along(nodes)) {
          tq <- xk * nodes[q]
          Hq <- build_h(tq, Xprev, spec_h)
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
      opt <- optim(theta0, fn = fn_k, method = "L-BFGS-B", control = list(maxit = maxit))
      a_hat <- if (m_g > 0) opt$par[seq_len(m_g)] else numeric(0)
      b_hat <- opt$par[(m_g + 1L):length(opt$par)]
      coeffs[[k]] <- list(alpha = a_hat, beta = b_hat)
    }
    S <- list(algo = "crossterm", mu = mu, sigma = sigma, coeffs = coeffs,
              deg_g = as.integer(deg_g), spec_h = spec_h,
              gl_nodes = nodes, gl_weights = weights, Hmax = Hmax)
    class(S) <- "ttm_cross_term"
  })[["elapsed"]]

  time_pred <- system.time({ invisible(predict_ttm(S, X_te, type = "logdensity_by_dim")) })[["elapsed"]]
  list(S = S,
       NLL_train = mean(-predict_ttm(S, X_tr, type = "logdensity")),
       NLL_test  = mean(-predict_ttm(S, X_te,  type = "logdensity")),
       time_train = time_train, time_pred = time_pred)
}

# S3 predict wrapper for cross-term TTM
predict.ttm_cross_term <- function(object, newdata, type = c("logdensity_by_dim", "logdensity"), ...) {
  predict_ttm(object, newdata, match.arg(type))
}

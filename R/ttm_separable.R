# Separable triangular map: S_k(x) = g_k(x_<k>) + f_k(x_k)
# Maps-from-samples (Eq. 38/39). Uses unified bases and core predictor.

source(file.path("R", "ttm_bases.R"))
source(file.path("R", "ttm_core.R"))

.std_stats <- function(X) {
  mu <- colMeans(X)
  sigma <- apply(X, 2, sd) + .Machine$double.eps
  list(mu = mu, sigma = sigma)
}

.standardize <- function(X, mu, sigma) {
  sweep(sweep(X, 2, mu, "-"), 2, sigma, "/")
}

fit_ttm_separable <- function(data, degree_g = 2L, lambda = 0.0, eps = 1e-6, seed = 42) {
  set.seed(seed)
  S_in <- if (is.list(data) && !is.null(data$X_tr) && !is.null(data$X_te)) data else {
    stopifnot(is.matrix(data))
    if (!exists("split_data")) source("02_split.R")
    list(X_tr = split_data(data, seed)$X_tr, X_te = split_data(data, seed)$X_te)
  }
  X_tr <- S_in$X_tr; X_te <- S_in$X_te
  K <- ncol(X_tr); N <- nrow(X_tr)

  time_train <- system.time({
    st <- .std_stats(X_tr); mu <- st$mu; sigma <- st$sigma
    Xs <- .standardize(X_tr, mu, sigma)
    coeffs <- vector("list", K)
    I_N <- diag(N)
    for (k in seq_len(K)) {
      x_prev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
      xk <- Xs[, k]
      Phi_non <- if (ncol(x_prev) > 0) build_g(x_prev, deg = degree_g) else if (degree_g == 0L) matrix(1, N, 1L) else matrix(0, N, 0)
      # For deg_g==0, use linear-only f and closed-form c for exact marginal equivalence
      if (degree_g == 0L) {
        # Exact marginal equivalence: use same closed-form as marginal per-dim
        v <- stats::var(xk) + 1e-12
        c_mon <- max(eps, 1 / sqrt(v))
        c_non <- matrix(-mean(xk) * c_mon, ncol = 1L)
        coeffs[[k]] <- list(c_non = c_non, c_mon = c_mon)
        next
      }
      Phi_mon <- build_f(xk)
      B <- d_build_f(xk)
      m_non <- ncol(Phi_non); m_mon <- ncol(Phi_mon)
      if (m_non > 0) {
        M <- solve(crossprod(Phi_non) + lambda * diag(m_non), t(Phi_non))
        A <- (I_N - Phi_non %*% M) %*% Phi_mon
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
      coeffs[[k]] <- list(c_non = c_non, c_mon = c_mon)
    }
    S <- list(algo = "separable", mu = mu, sigma = sigma, coeffs = coeffs, degree_g = degree_g)
    class(S) <- "ttm_separable"
    model <<- S
  })[["elapsed"]]

  time_pred <- system.time({ invisible(predict_ttm(model, X_te, type = "logdensity_by_dim")) })[["elapsed"]]

  list(S = model,
       NLL_train = mean(-predict_ttm(model, X_tr, type = "logdensity")),
       NLL_test  = mean(-predict_ttm(model, X_te,  type = "logdensity")),
       time_train = time_train, time_pred = time_pred)
}

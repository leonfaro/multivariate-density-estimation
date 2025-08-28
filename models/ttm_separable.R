# Triangular Transport Map - Separable Modul
# Basis-R Implementierung nach Appendix A

# Hilfsfunktionen ------------------------------------------------------------

if (!exists(".standardizeData")) {
  .standardizeData <- function(X) {
    mu <- colMeans(X)
    sigma <- apply(X, 2, sd) + .Machine$double.eps
    X_tilde <- sweep(X, 2, mu, "-")
    X_tilde <- sweep(X_tilde, 2, sigma, "/")
    list(X = X_tilde, mu = mu, sigma = sigma)
  }
}

if (!exists(".standardize")) {
  .standardize <- function(S, X) {
    X <- sweep(X, 2, S$mu, "-")
    sweep(X, 2, S$sigma, "/")
  }
}

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

basis_f <- function(x) cbind(x, erf(x))

dbasis_f <- function(x) {
  cbind(rep(1, length(x)), 2 / sqrt(pi) * exp(-x^2))
}

basis_g <- function(X, deg) {
  if (ncol(X) == 0L) {
    return(matrix(0, nrow = nrow(X), ncol = 0))
  }
  N <- nrow(X)
  out <- matrix(1, N, 1)
  for (j in seq_len(ncol(X))) {
    xj <- X[, j]
    for (d in seq_len(deg)) {
      out <- cbind(out, xj^d)
    }
  }
  out
}

# Exportierte Funktionen -----------------------------------------------------

trainSeparableMap <- function(X_or_path, degree_g = 2, lambda = 1e-3, eps = 1e-6, seed = 42) {
  set.seed(seed)
  S_in <- if (is.character(X_or_path)) readRDS(X_or_path) else X_or_path
  stopifnot(is.list(S_in))
  X_tr <- S_in$X_tr
  X_te  <- S_in$X_te

  time_train <- system.time({
    std <- .standardizeData(X_tr)
    X_tr_std <- std$X
    mu <- std$mu
    sigma <- std$sigma
    K <- ncol(X_tr_std)
    coeffs <- vector("list", K)
    N <- nrow(X_tr_std)
    for (k in seq_len(K)) {
      x_prev <- if (k > 1) X_tr_std[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
      xk <- X_tr_std[, k]
      P_non <- if (k > 1) basis_g(x_prev, degree_g) else matrix(0, N, 0)
      P_mon <- basis_f(xk)
      B <- dbasis_f(xk)
      stopifnot(nrow(P_mon) == N, nrow(P_non) == N, nrow(B) == N)
      m_non <- ncol(P_non)
      m_mon <- ncol(P_mon)
      if (m_non > 0) {
        M <- solve(crossprod(P_non) + lambda * diag(m_non), t(P_non))
        A <- (diag(N) - P_non %*% M) %*% P_mon
        D <- M %*% P_mon
      } else {
        M <- matrix(0, 0, N)
        A <- P_mon
        D <- matrix(0, 0, m_mon)
      }
      fn <- function(c) {
        r <- A %*% c
        Bc <- B %*% c
        if (any(Bc <= 0)) return(Inf)
        q <- D %*% c
        0.5 * sum(r^2) - sum(log(Bc)) + (lambda / 2) * (sum(q^2) + sum(c^2))
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
      c_non <- if (m_non > 0) -M %*% (P_mon %*% c_mon) else numeric(0)
      coeffs[[k]] <- list(c_non = c_non, c_mon = c_mon)
    }
    S_map <- list(
      mu = mu,
      sigma = sigma,
      coeffs = coeffs,
      degree_g = degree_g,
      order = seq_len(K)
    )
    class(S_map) <- "ttm_separable"
  })[["elapsed"]]
  time_pred <- system.time({
    predict(S_map, X_te, "logdensity_by_dim")
  })[["elapsed"]]

  list(
    S = S_map,
    NLL_train = NLL_set(S_map, X_tr),
    NLL_test = NLL_set(S_map, X_te),
    stderr_test = SE_set(S_map, X_te),
    time_train = time_train,
    time_pred = time_pred
  )
}

predict.ttm_separable <- function(object, newdata,
                                  type = c("logdensity_by_dim", "logdensity")) {
  type <- tryCatch(match.arg(type), error = function(e) stop("unknown type"))
  Xs <- .standardize(object, newdata)
  N <- nrow(Xs)
  K <- ncol(Xs)
  Z <- matrix(0, N, K)
  LJ <- matrix(0, N, K)
  for (k in seq_len(K)) {
    x_prev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
    xk <- Xs[, k]
    P_non <- if (k > 1) basis_g(x_prev, object$degree_g) else matrix(0, N, 0)
    P_mon <- basis_f(xk)
    B <- dbasis_f(xk)
    c_non <- object$coeffs[[k]]$c_non
    c_mon <- object$coeffs[[k]]$c_mon
    gk <- if (ncol(P_non) > 0) as.numeric(P_non %*% c_non) else rep(0, N)
    fk <- as.numeric(P_mon %*% c_mon)
    Z[, k] <- gk + fk
    deriv <- (B %*% c_mon) / object$sigma[k]
    LJ[, k] <- log(as.numeric(deriv))
  }
  C <- -0.5 * log(2 * pi)
  LD <- (-0.5) * (Z^2) + C + LJ
  # Invariants and numeric sanity
  stopifnot(is.matrix(LD), nrow(LD) == N, ncol(LD) == K)
  if (!all(is.finite(LD))) stop("Non-finite values in separable logdensity_by_dim")
  if (type == "logdensity_by_dim") {
    LD_joint <- rowSums(LD)
    if (!all(is.finite(LD_joint))) stop("Non-finite joint logdensity (separable)")
    stopifnot(max(abs(LD_joint - rowSums(LD))) <= 1e-10)
    LD
  } else {
    LD_joint <- rowSums(LD)
    if (!all(is.finite(LD_joint))) stop("Non-finite joint logdensity (separable)")
    stopifnot(max(abs(LD_joint - rowSums(LD))) <= 1e-10)
    LD_joint
  }
}

NLL_set <- function(S, X) {
  mean(-rowSums(predict(S, X, "logdensity_by_dim")))
}

SE_set <- function(S, X) {
  v <- rowSums(-predict(S, X, "logdensity_by_dim"))
  stats::sd(v) / sqrt(length(v))
}

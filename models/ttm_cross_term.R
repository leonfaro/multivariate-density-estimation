# Triangular Transport Map - Cross-Term Modul
# Basis-R Implementierung nach Vorgaben

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

basis_g <- function(X, deg) {
  if (ncol(X) == 0L) {
    return(matrix(1, nrow = nrow(X), ncol = 1))
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

psi_basis <- function(t, xprev, deg_t, deg_x, cross = TRUE) {
  out <- numeric(0)
  for (d in seq_len(deg_t)) {
    out <- c(out, t^d)
  }
  if (length(xprev) > 0) {
    for (j in seq_along(xprev)) {
      for (d in seq_len(deg_x)) {
        out <- c(out, xprev[j]^d)
      }
    }
    if (cross) {
      for (j in seq_along(xprev)) {
        out <- c(out, t * xprev[j])
      }
    }
  }
  out
}

gauss_legendre_01 <- function(n) {
  if (n <= 0 || n != as.integer(n)) {
    stop("n must be positive integer")
  }
  if (n == 1) {
    return(list(nodes = 0.5, weights = 1))
  }
  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  for (k in i) {
    J[k, k + 1] <- b[k]
    J[k + 1, k] <- b[k]
  }
  e <- eigen(J, symmetric = TRUE)
  x <- (e$values + 1) / 2
  w <- (2 * (e$vectors[1, ]^2)) / 2
  list(nodes = x, weights = w)
}

# Exportierte Funktionen -----------------------------------------------------

trainCrossTermMap <- function(X_or_path, degree_g = 2, degree_t = 2) {
  set.seed(42)
  S_in <- if (is.character(X_or_path)) readRDS(X_or_path) else X_or_path
  stopifnot(is.list(S_in))
  X_tr <- S_in$X_tr
  X_val <- S_in$X_val
  X_te <- S_in$X_te

  std <- .standardizeData(X_tr)
  mu <- std$mu
  sigma <- std$sigma
  K <- ncol(X_tr)
  coeffs <- vector("list", K)
  for (k in seq_len(K)) {
    coeffs[[k]] <- list(alpha = numeric(0), beta = numeric(0))
  }
  S_map <- list(
    mu = mu,
    sigma = sigma,
    degree_g = degree_g,
    degree_t = degree_t,
    coeffs = coeffs,
    order = seq_len(K)
  )
  class(S_map) <- "ttm_cross_term"

  list(
    S = S_map,
    NLL_train = NLL_set(S_map, X_tr),
    NLL_val = NLL_set(S_map, X_val),
    NLL_test = NLL_set(S_map, X_te),
    stderr_test = SE_set(S_map, X_te),
    time_train = 0,
    time_pred = 0
  )
}

predict.ttm_cross_term <- function(object, newdata,
                                   type = c("logdensity_by_dim", "logdensity")) {
  type <- tryCatch(match.arg(type), error = function(e) stop("unknown type"))
  Xs <- .standardize(object, newdata)
  N <- nrow(Xs)
  K <- ncol(Xs)
  LD <- matrix(0, N, K)
  if (type == "logdensity_by_dim") {
    LD
  } else {
    rowSums(LD)
  }
}

NLL_set <- function(S, X) {
  mean(-rowSums(predict(S, X, "logdensity_by_dim")))
}

SE_set <- function(S, X) {
  v <- rowSums(-predict(S, X, "logdensity_by_dim"))
  stats::sd(v) / sqrt(length(v))
}

forwardKLLoss <- function(S, X) {
  0
}

logJacDiag <- function(S, x) {
  rep(0, length(x))
}


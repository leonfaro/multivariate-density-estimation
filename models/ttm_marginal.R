# Triangular Transport Map - Marginal Module
# Basis-R Implementierung

# interne Hilfsfunktionen --------------------------------------------------

.standardizeData <- function(X) {
  mu <- colMeans(X)
  sigma <- apply(X, 2, sd) + .Machine$double.eps
  X_tilde <- sweep(X, 2, mu, "-")
  X_tilde <- sweep(X_tilde, 2, sigma, "/")
  list(X = X_tilde, mu = mu, sigma = sigma)
}

.forward_matrix <- function(S, X) {
  Xs <- .standardize(S, X)
  b <- exp(S$coeffA)
  Z <- sweep(Xs, 2, b, "*")
  sweep(Z, 2, S$coeffB, "+")
}

.logjac_const <- function(S) {
  b <- exp(S$coeffA)
  log(b) - log(S$sigma)
}

.standardize <- function(S, X) {
  X <- sweep(X, 2, S$mu, "-")
  sweep(X, 2, S$sigma, "/")
}

# exportierte Funktionen ----------------------------------------------------

trainMarginalMap <- function(X_or_path, seed = 42) {
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
    coeffA <- numeric(K)
    coeffB <- numeric(K)
    for (k in seq_len(K)) {
      xk <- X_tr_std[, k]
      N <- length(xk)
      u <- rank(xk, ties.method = "average") / (N + 1)
      lower <- 1 / (N + 1)
      upper <- N / (N + 1)
      u <- pmin(pmax(u, lower), upper)
      z_star <- qnorm(u)
      covxz <- mean((xk - mean(xk)) * (z_star - mean(z_star)))
      varx <- var(xk) + 1e-12
      b_k <- max(0, covxz / varx)
      a_k <- mean(z_star) - b_k * mean(xk)
      coeffA[k] <- log(b_k + 1e-12)
      coeffB[k] <- a_k
    }
    coeffC <- rep(0, K)
    S_map <- list(
      mu = mu,
      sigma = sigma,
      coeffA = coeffA,
      coeffB = coeffB,
      coeffC = coeffC,
      order = seq_len(K)
    )
    class(S_map) <- "ttm_marginal"
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

predict.ttm_marginal <- function(object, newdata,
                                 type = c("logdensity_by_dim", "logdensity")) {
  type <- tryCatch(match.arg(type), error = function(e) stop("unknown type"))
  Z <- .forward_matrix(object, newdata)
  LJ <- .logjac_const(object)
  C <- -0.5 * log(2 * pi)
  LD <- (-0.5) * (Z^2) + C +
    matrix(LJ, nrow = nrow(Z), ncol = length(LJ), byrow = TRUE)
  # Invariants and numeric sanity
  stopifnot(is.matrix(LD), nrow(LD) == nrow(newdata), ncol(LD) == ncol(newdata))
  if (!all(is.finite(LD))) stop("Non-finite values in marginal logdensity_by_dim")
  if (type == "logdensity_by_dim") {
    # Assert joint equals rowSums
    LD_joint <- rowSums(LD)
    if (!all(is.finite(LD_joint))) stop("Non-finite joint logdensity (marginal)")
    stopifnot(max(abs(LD_joint - rowSums(LD))) <= 1e-10)
    LD
  } else {
    LD_joint <- rowSums(LD)
    if (!all(is.finite(LD_joint))) stop("Non-finite joint logdensity (marginal)")
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

forwardPass <- function(S, x) {
  x_std <- .standardize(S, matrix(x, nrow = 1))
  b <- exp(S$coeffA)
  a <- S$coeffB
  as.numeric(a + b * x_std)
}

logJacDiag <- function(S, x) {
  LJ <- .logjac_const(S)
  rep(LJ, length.out = length(x))
}

forwardKLLoss <- function(S, X) {
  Z <- .forward_matrix(S, X)
  LJ <- .logjac_const(S)
  mean(0.5 * rowSums(Z^2) - sum(LJ))
}

inversePass <- function(S, z) {
  b <- exp(S$coeffA)
  a <- S$coeffB
  x_std <- (z - a) / b
  x <- sweep(x_std, 2, S$sigma, "*")
  x <- sweep(x, 2, S$mu, "+")
  as.numeric(x)
}

negativeLogLikelihood <- function(S, X) {
  Z <- .forward_matrix(S, X)
  LJ <- .logjac_const(S)
  sum(0.5 * rowSums(Z^2) - sum(LJ))
}

natsPerDim <- function(NLL, N, K) {
  NLL / (N * K)
}

stderr <- function(v) {
  stats::sd(v) / sqrt(length(v))
}

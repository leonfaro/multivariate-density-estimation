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

.computeRowwiseLosses <- function(S, X_set) {
  b <- exp(S$coeffA)
  a <- S$coeffB
  z <- sweep(X_set, 2, b, "*") +
    matrix(a, nrow = nrow(X_set), ncol = length(a), byrow = TRUE)
  0.5 * rowSums(z^2) - sum(S$coeffA)
}

.standardize <- function(S, X) {
  X <- sweep(X, 2, S$mu, "-")
  sweep(X, 2, S$sigma, "/")
}

# exportierte Funktionen ----------------------------------------------------

trainMarginalMap <- function(X_or_path) {
  set.seed(42)
  S_in <- if (is.character(X_or_path)) readRDS(X_or_path) else X_or_path
  stopifnot(is.list(S_in))
  X_tr <- S_in$X_tr
  X_val <- S_in$X_val
  X_te  <- S_in$X_te

  time_train <- system.time({
    std <- .standardizeData(X_tr)
    X_tr_std <- std$X
    mu <- std$mu
    sigma <- std$sigma
    X_val_std <- sweep(sweep(X_val, 2, mu, "-"), 2, sigma, "/")
    X_te_std  <- sweep(sweep(X_te,  2, mu, "-"), 2, sigma, "/")
    K <- ncol(X_tr_std)
    coeffA <- numeric(K)
    coeffB <- numeric(K)
    for (k in seq_len(K)) {
      xk <- X_tr_std[, k]
      u <- rank(xk, ties.method = "average") / (length(xk) + 1)
      lower <- 1 / (length(xk) + 1)
      upper <- length(xk) / (length(xk) + 1)
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
    loss_train_vec <- .computeRowwiseLosses(S_map, X_tr_std)
    loss_val_vec <- .computeRowwiseLosses(S_map, X_val_std)
  })[["elapsed"]]

  loss_test_time <- system.time({
    loss_test_vec <- .computeRowwiseLosses(S_map, X_te_std)
  })[["elapsed"]]

  list(
    S = S_map,
    NLL_train = mean(loss_train_vec),
    NLL_val = mean(loss_val_vec),
    NLL_test = mean(loss_test_vec),
    stderr_test = stderr(loss_test_vec),
    time_train = time_train,
    time_pred = loss_test_time
  )
}

predict.ttm_marginal <- function(object, newdata,
                                 type = c("logdensity_by_dim", "logdensity")) {
  type <- match.arg(type)
  X_std <- .standardize(object, newdata)
  b <- exp(object$coeffA)
  a <- object$coeffB
  z <- sweep(X_std, 2, b, "*") +
    matrix(a, nrow = nrow(X_std), ncol = length(a), byrow = TRUE)
  log_diag <- matrix(object$coeffA,
                     nrow = nrow(X_std), ncol = length(a), byrow = TRUE)
  log_dens_dim <- -0.5 * z^2 + log_diag
  if (type == "logdensity_by_dim") {
    log_dens_dim
  } else {
    rowSums(log_dens_dim)
  }
}

forwardPass <- function(S, x) {
  x_std <- .standardize(S, matrix(x, nrow = 1))
  b <- exp(S$coeffA)
  a <- S$coeffB
  as.numeric(a + b * x_std)
}

logJacDiag <- function(S, x) {
  rep(S$coeffA, length.out = length(x))
}

forwardKLLoss <- function(S, X) {
  X_std <- .standardize(S, X)
  mean(.computeRowwiseLosses(S, X_std))
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
  X_std <- .standardize(S, X)
  sum(.computeRowwiseLosses(S, X_std))
}

natsPerDim <- function(NLL, N, K) {
  NLL / (N * K)
}

stderr <- function(v) {
  stats::sd(v) / sqrt(length(v))
}

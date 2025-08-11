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

.shuffleOrdering <- function(K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(seq_len(K))
}

.updateCoeffsMarginal <- function(S, X_batch) {
  d <- length(S$order)
  for (j in seq_len(d)) {
    k <- S$order[j]
    xk <- X_batch[, k]
    u <- rank(xk, ties.method = "average") / (length(xk) + 1)
    z_star <- qnorm(u)
    covxz <- mean((xk - mean(xk)) * (z_star - mean(z_star)))
    varx <- var(xk) + 1e-12
    b_star <- max(0, covxz / varx)
    a_star <- mean(z_star) - b_star * mean(xk)
    S$coeffA[k] <- log(b_star + 1e-12)
    S$coeffB[k] <- a_star
  }
  S
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

  X_all <- rbind(X_tr, X_val, X_te)
  std <- .standardizeData(X_all)
  X_std <- std$X
  n_tr <- nrow(X_tr)
  n_val <- nrow(X_val)
  X_train <- X_std[seq_len(n_tr), , drop = FALSE]
  X_val_std <- X_std[seq_len(n_val) + n_tr, , drop = FALSE]
  X_test <- X_std[(n_tr + n_val + 1):nrow(X_std), , drop = FALSE]

  K <- ncol(X_train)
  order <- .shuffleOrdering(K)
  S_map <- list(
    mu = std$mu,
    sigma = std$sigma,
    coeffA = rep(0, K),
    coeffB = rep(0, K),
    coeffC = rep(0, K),
    order = order
  )
  class(S_map) <- "ttm_marginal"

  best_val <- Inf
  best_state <- S_map
  best_epoch <- 0L
  best_train <- Inf
  patience <- 0L
  T_max <- 100L
  P <- 10L

  time_train <- system.time({
    for (epoch in seq_len(T_max)) {
      S_map <- .updateCoeffsMarginal(S_map, X_train)
      NLL_train <- mean(.computeRowwiseLosses(S_map, X_train))
      NLL_val <- mean(.computeRowwiseLosses(S_map, X_val_std))

      if (NLL_val < best_val - 1e-6) {
        best_val <- NLL_val
        best_state <- S_map
        best_epoch <- epoch
        best_train <- NLL_train
        patience <- 0L
      } else {
        patience <- patience + 1L
      }
      if (patience > P) break
    }
  })[["elapsed"]]

  S_map <- best_state

  time_pred <- system.time({
    loss_test_vec <- .computeRowwiseLosses(S_map, X_test)
  })[["elapsed"]]

  NLL_test <- mean(loss_test_vec)
  stderr_test <- stderr(loss_test_vec)

  list(
    S = S_map,
    NLL_train = best_train,
    NLL_val = best_val,
    NLL_test = NLL_test,
    stderr_test = stderr_test,
    time_train = time_train,
    time_pred = time_pred
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


# Triangular Transport Maps - Generic Function Skeleton

# --- 1 Basis-Datendienst -------------------------------------------------------

standardizeData <- function(X) {
  mu <- colMeans(X)
  sigma <- apply(X, 2, sd) + .Machine$double.eps
  X_tilde <- sweep(X, 2, mu, '-')
  X_tilde <- sweep(X_tilde, 2, sigma, '/')
  list(X = X_tilde, mu = mu, sigma = sigma)
}

sampleReference <- function(N, d) {
  matrix(rnorm(N * d), nrow = N, ncol = d)
}

shuffleOrdering <- function(K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(seq_len(K))
}

# --- 2 Map-Abstraktion ---------------------------------------------------------

MapStruct <- function(type, coeffA, coeffB, coeffC, basisF, basisG, basisH) {
  structure(
    list(
      type = type,
      coeffA = coeffA,
      coeffB = coeffB,
      coeffC = coeffC,
      basisF = basisF,
      basisG = basisG,
      basisH = basisH
    ),
    class = 'MapStruct'
  )
}

# --- 3 Elementare Algebra-Funktionen ------------------------------------------

evaluateMap <- function(S, x) {
  forwardPass(S, x)
}

jacDiag <- function(S, x) {
  d <- length(x)
  diag <- numeric(d)
  for (k in seq_len(d)) {
    if (S$type == 'marginal') {
      diag[k] <- S$basisF[[k]]$deriv(x[k], S$coeffA[[k]])
    } else if (S$type == 'separable') {
      diag[k] <- S$basisF[[k]]$deriv(x[k], S$coeffA[[k]])
    } else {
      diag[k] <- exp(S$basisH[[k]](x[k], x[seq_len(k - 1)], S$coeffC[[k]]))
    }
  }
  diag
}

detJacobian <- function(diag) {
  prod(diag)
}

logDetJacobian <- function(diag) {
  sum(log(diag))
}

forwardKLLoss <- function(S, X) {
  total <- 0
  for (i in seq_len(nrow(X))) {
    z <- forwardPass(S, X[i, ])
    diag <- jacDiag(S, X[i, ])
    total <- total + 0.5 * sum(z^2) - logDetJacobian(diag)
  }
  total / nrow(X)
}

# --- 4 Numerische Grundbausteine ----------------------------------------------

basisEval1D <- function(basisSet, x) {
  basisSet(x)
}

basisEvalKD <- function(basisSet, vec) {
  do.call(basisSet, as.list(vec))
}

monotoneIntegrator <- function(h, t0, t) {
  integrate(function(s) exp(h(s)), lower = t0, upper = t)$value
}

rootFind1D <- function(fun, target) {
  uniroot(function(t) fun(t) - target, interval = c(-10, 10))$root
}

optimStep <- function(params, grad, lr) {
  params - lr * grad
}

batchIterator <- function(X, B) {
  N <- nrow(X)
  idx <- split(sample(seq_len(N)), ceiling(seq_along(seq_len(N)) / B))
  lapply(idx, function(i) X[i, , drop = FALSE])
}

# --- 5 Hilfs-Workflows --------------------------------------------------------

forwardPass <- function(S, x) {
  d <- length(x)
  z <- numeric(d)
  for (k in seq_len(d)) {
    if (S$type == 'marginal') {
      z[k] <- S$basisF[[k]](x[k], S$coeffA[[k]])
    } else if (S$type == 'separable') {
      z[k] <- S$basisG[[k]](x[seq_len(k - 1)], S$coeffB[[k]]) +
        S$basisF[[k]](x[k], S$coeffA[[k]])
    } else {
      z[k] <- S$basisG[[k]](x[seq_len(k - 1)], S$coeffB[[k]]) +
        monotoneIntegrator(
          function(s) S$basisH[[k]](s, x[seq_len(k - 1)], S$coeffC[[k]]),
          0,
          x[k]
        )
    }
  }
  z
}

inversePass <- function(S, z, tol = 1e-8) {
  d <- length(z)
  x <- numeric(d)
  for (k in seq_len(d)) {
    fun <- function(t) {
      temp <- x
      temp[k] <- t
      forwardPass(S, temp)[k]
    }
    x[k] <- rootFind1D(fun, z[k])
  }
  x
}

lossFull <- function(S, X) {
  forwardKLLoss(S, X)
}

# --- 6 Generische Train-Routine ----------------------------------------------

trainMap <- function(S, X, nIter, lr, batchSize) {
  for (iter in seq_len(nIter)) {
    batches <- batchIterator(X, batchSize)
    for (B in batches) {
      grad <- lapply(S$coeffA, function(a) rep(0, length(a)))
      S$coeffA <- mapply(optimStep, S$coeffA, grad, MoreArgs = list(lr = lr), SIMPLIFY = FALSE)
    }
  }
  S
}

# --- 7 Evaluations-Utilities --------------------------------------------------

negativeLogLikelihood <- function(S, Xtest) {
  apply(Xtest, 1, function(x) {
    z <- forwardPass(S, x)
    diag <- jacDiag(S, x)
    0.5 * sum(z^2) - logDetJacobian(diag)
  })
}

natsPerDim <- function(NLL, N) {
  sum(NLL) / (length(NLL) * N)
}

stderr <- function(values) {
  stats::sd(values) / sqrt(length(values))
}


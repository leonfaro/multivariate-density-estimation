# TTM - Marginal Map Trainer
# arbeitet komplett in Basis-R
source("models/ttm_base.R")

## Reproduzierbarkeit -------------------------------------------------------
set.seed(42)

# Globale Hyperparameter ---------------------------------------------------
lr0 <- 0.01
T_max <- 100L
P <- 10L
decay <- 1.0

# Datenladefunktion -------------------------------------------------------
loadCSV <- function(filepath) {
  as.matrix(read.csv(filepath))
}

## Hilfsfunktionen ----------------------------------------------------------

initializeCoeffs <- function(S) {
  d <- length(S$order)
  S$coeffA <- vector("list", d)
  S$coeffB <- vector("list", d)
  S$coeffC <- vector("list", d)
  for (k in seq_len(d)) {
    # log(1.0) = 0   ⇒ f' = exp(0) = 1 > 0
    S$coeffA[[k]] <- log(1.0)
    S$coeffB[[k]] <- 0
    S$coeffC[[k]] <- 0
  }
  S
}

updateCoeffsMarginal <- function(S, X_batch, lr) {
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

    S$coeffA[[k]] <- log(b_star + 1e-12)
    S$coeffB[[k]] <- a_star
    S$coeffC[[k]] <- 0
  }
  S
}

computeRowwiseLosses <- function(S, X_set) {
  losses <- numeric(nrow(X_set))
  for (i in seq_len(nrow(X_set))) {
    z <- forwardPass(S, X_set[i, ])
    ell <- logJacDiag(S, X_set[i, ])
    losses[i] <- 0.5 * sum(z^2) - sum(ell)
  }
  losses
}

## Hauptfunktion ------------------------------------------------------------

trainMarginalMap <- function(filepath) {
  X_raw <- loadCSV(filepath)
  std_res <- standardizeData(X_raw)
  X_std  <- std_res[[1]]
  N <- nrow(X_std)
  d <- ncol(X_std)

  idx <- shuffleOrdering(N)
  train_idx <- idx[1:floor(0.8 * N)]
  val_idx   <- idx[(floor(0.8 * N) + 1):floor(0.9 * N)]
  test_idx  <- idx[(floor(0.9 * N) + 1):N]

  X_train <- X_std[train_idx, , drop = FALSE]
  X_val   <- X_std[val_idx, , drop = FALSE]
  X_test  <- X_std[test_idx, , drop = FALSE]

  S <- MapStruct(type = "marginal")
  S <- setOrdering(S, shuffleOrdering(d))
  S <- initializeCoeffs(S)
  if (is.null(S$basisF)) {
    S$basisF <- vector("list", d)
  }
  for (k in seq_len(d)) {
    S$basisF[[k]] <- function(x) x
  }

  best_val <- Inf
  best_state <- S
  best_epoch <- 0L
  patience <- 0L
  lr <- lr0

  for (epoch in seq_len(T_max)) {
    S <- updateCoeffsMarginal(S, X_train, lr)
    NLL_train <- mean(computeRowwiseLosses(S, X_train))
    NLL_val <- mean(computeRowwiseLosses(S, X_val))

    if (NLL_val < best_val - 1e-6) {
      best_val <- NLL_val
      best_state <- S
      best_epoch <- epoch
      patience <- 0L
    } else {
      patience <- patience + 1L
    }
    if (patience > P) break
    lr <- lr * decay
    if (epoch %% 10 == 0) {
      message(epoch, ": val NLL = ", round(NLL_val, 4))
    }
  }

  S <- best_state
  loss_test_vec <- computeRowwiseLosses(S, X_test)
  NLL_test <- mean(loss_test_vec)
  stderr_test <- stderr(loss_test_vec)

  list(
    S = S,
    best_epoch = best_epoch,
    NLL_train = NLL_train,
    NLL_val = best_val,
    NLL_test = NLL_test,
    stderr_test = stderr_test
  )
}


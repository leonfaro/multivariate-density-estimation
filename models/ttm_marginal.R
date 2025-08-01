# TTM - Marginal Map Trainer
# arbeitet komplett in Basis-R

## Reproduzierbarkeit -------------------------------------------------------
set.seed(42)

# Globale Hyperparameter ---------------------------------------------------
lr0 <- 0.01
T_max <- 100L
P <- 10L
decay <- 1.0

## Hilfsfunktionen ----------------------------------------------------------

linearBasis <- function(S, idx) {
  f <- function(x, theta) {
    S$coeffB[[idx]] + exp(theta) * x
  }
  attr(f, "deriv") <- function(x, theta) rep(exp(theta), length(x))
  f
}

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

trainMarginalMap <- function(S) {
  stopifnot(is.list(S))
  X_tr  <- S$X_tr
  X_val <- S$X_val
  X_te  <- S$X_te

  X_all <- rbind(X_tr, X_val, X_te)
  std_res <- standardizeData(X_all)
  X_std <- std_res$X
  n_tr <- nrow(X_tr)
  n_val <- nrow(X_val)
  X_train <- X_std[seq_len(n_tr), , drop = FALSE]
  X_val   <- X_std[seq_len(n_val) + n_tr, , drop = FALSE]
  X_test  <- X_std[(n_tr + n_val + 1):nrow(X_std), , drop = FALSE]

  d <- ncol(X_train)
  S_map <- MapStruct(type = "marginal")
  S_map <- setOrdering(S_map, shuffleOrdering(d))
  S_map <- initializeCoeffs(S_map)
  S_map$basisF <- vector("list", d)
  for (k in seq_len(d)) {
    S_map$basisF[[k]] <- linearBasis(S_map, k)
  }

  best_val <- Inf
  best_state <- S_map
  best_epoch <- 0L
  best_train <- Inf
  patience <- 0L
  lr <- lr0

  for (epoch in seq_len(T_max)) {
    S_map <- updateCoeffsMarginal(S_map, X_train, lr)
    NLL_train <- mean(computeRowwiseLosses(S_map, X_train))
    NLL_val <- mean(computeRowwiseLosses(S_map, X_val))

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
    lr <- lr * decay
    if (epoch %% 10 == 0) {
      message(epoch, ": val NLL = ", round(NLL_val, 4))
    }
  }

  S_map <- best_state
  loss_test_vec <- computeRowwiseLosses(S_map, X_test)
  NLL_test <- mean(loss_test_vec)
  stderr_test <- stderr(loss_test_vec)

  list(
    S = S_map,
    best_epoch = best_epoch,
    NLL_train = best_train,
    NLL_val = best_val,
    NLL_test = NLL_test,
    stderr_test = stderr_test
  )
}

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
  stop("updateCoeffsMarginal() noch nicht implementiert")
}

computeRowwiseLosses <- function(S, X_set) {
  stop("computeRowwiseLosses() noch nicht implementiert")
}

## Hauptfunktion ------------------------------------------------------------

trainMarginalMap <- function(filepath) {
  X_raw <- loadCSV(filepath)
  std_res <- standardizeData(X_raw)
  X_std  <- std_res[[1]]
  mu     <- std_res[[2]]
  sigma  <- std_res[[3]]
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

  lr <- lr0
  Tmax <- T_max
  Patience <- P
  lr_decay <- decay

  list(S = S,
       X_train = X_train,
       X_val = X_val,
       X_test = X_test,
       lr = lr,
       Tmax = Tmax,
       P = Patience,
       decay = lr_decay)
}


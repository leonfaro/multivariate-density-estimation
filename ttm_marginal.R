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

initializeCoeffs <- function(S, d) {
  S$coeffA <- lapply(seq_len(d), function(k) 0)
  S$coeffB <- vector("list", d)
  S$coeffC <- vector("list", d)
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
  X_std  <- std_res$X
  mu     <- std_res$mu
  sigma  <- std_res$sigma
  N <- nrow(X_std)
  d <- ncol(X_std)

  idx <- shuffleOrdering(N)
  train_idx <- idx[1:floor(0.8 * N)]
  val_idx   <- idx[(floor(0.8 * N) + 1):floor(0.9 * N)]
  test_idx  <- idx[(floor(0.9 * N) + 1):N]

  X_train <- X_std[train_idx, , drop = FALSE]
  X_val   <- X_std[val_idx, , drop = FALSE]
  X_test  <- X_std[test_idx, , drop = FALSE]

  S <- MapStruct(type = "marginal",
                 coeffA = NULL,
                 coeffB = NULL,
                 coeffC = NULL,
                 basisF = vector("list", d),
                 basisG = NULL,
                 basisH = NULL)

  S$order <- shuffleOrdering(d)
  S <- initializeCoeffs(S, d)

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


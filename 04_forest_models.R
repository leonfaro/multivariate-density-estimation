
## -------------------------------------------------------------------
## Overview of this script
##
## Goal: Train a transformation forest estimator on observation
##       vectors \(x^{(i)} \in \mathbb R^K\) and evaluate the
##       sequential log-density on a separate test sample.  Both
##       training and testing rely exclusively on the matrices
##       `X_pi_train` and `X_pi_test` generated in
##       `02_generate_data.R`.
##
## - Output: transformation forest model `model` and matrix `LD_hat`
##   containing log-density contributions for the test data.
## - Algorithm:
##   1. For each coordinate \(k=1,\dots,K\) fit a marginal Box--Cox
##      model on the single column `X_pi_train[, k]`.
##   2. For every \(k>1\) fit a conditional transformation forest for
##      the response `X_pi_k` given all preceding coordinates
##      `X_pi_{<k}`.  The fit uses only these columns.
##   3. Predict log-density contributions on `X_pi_test` sequentially
##      and assemble them into `LD_hat`.
## -------------------------------------------------------------------
source("02_generate_data.R")

library(trtf)

## reproducibility for the forest fit
set.seed(2046)

fit_forest <- function(data) {
  ymod <- lapply(seq_len(ncol(data)), function(j)
    BoxCox(as.formula(paste0(names(data)[j], "~1")), data = data[, j, drop = FALSE])
  )
  forests <- vector("list", ncol(data) - 1L)
  for (k in 2:ncol(data)) {
    fm <- as.formula(paste0(
      names(data)[k], " ~ ",
      paste(names(data)[1:(k-1)], collapse = "+")
    ))
    forests[[k - 1L]] <- traforest(
      ymod[[k]], formula = fm, data = data[, seq_len(k)],
      ntree = 200,
      mtry = ceiling((k - 1) / 3),
      minbucket = 20,
      trace = TRUE
    )
  }
  structure(list(ymod = ymod, forests = forests), class = "mytrtf")
}

predict.mytrtf <- function(object, newdata,
                           type = c("logdensity", "distribution",
                                     "density", "trafo", "response")) {
  type <- match.arg(type)
  K <- length(object$ymod)
  ld <- matrix(NA_real_, nrow(newdata), K)
  ## k = 1: marginal model
  ld[, 1] <- predict(object$ymod[[1]], newdata = newdata[, 1, drop = FALSE],
                     type = type)
  for (k in 2:K) {
      frst <- object$forests[[k - 1L]]
      nd_sub <- newdata[, seq_len(k), drop = FALSE]
      resp <- variable.names(frst$model)[1]
      q <- nd_sub[[resp]]
      pr <- predict(frst, newdata = nd_sub, q = q,
                    type = type, simplify = FALSE)
      ld[, k] <- diag(do.call(cbind, pr))
  }
  stopifnot(all(is.finite(ld)))
  ld
}

model  <- fit_forest(as.data.frame(X_pi_train))
LD_hat <- predict(model, newdata = as.data.frame(X_pi_test), type = "logdensity")

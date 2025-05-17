# Fit transformation forests and evaluate joint log-likelihood
# This script trains conditional transformation forests for each
# distributional dimension and compares the resulting joint log-likelihood
# to the ground truth values computed earlier.
source("baseline_parametric.R")

library(trtf)
library(tram)

#' Fit a series of conditional transformation forests
#'
#' @param data Data frame of training observations.
#' @return Object with marginal transformation models and forests.
mytrtf <- function(data) {
  ymod <- lapply(names(data), function(y)
    BoxCox(as.formula(paste0(y, "~1")), data = data)
  )
  fm <- lapply(2:ncol(data), function(j)
    as.formula(paste0(
      names(data)[j], " ~ ",
      paste(names(data)[1:(j-1)], collapse = "+")
    ))
  )
  forests <- lapply(seq_along(fm), function(j)
    traforest(
      ymod[[j+1]], formula = fm[[j]], data = data,
      ntree = 200,
      mtry = ceiling((j-1)/3),
      minbucket = 20,
      trace = TRUE
    )
  )
  structure(list(ymod = ymod, forests = forests), class = "mytrtf")
}

 #' Predict log-density contributions from a fitted mytrtf object
 #'
 #' @param object Fitted model returned by `mytrtf`.
 #' @param newdata Data frame of observations to score.
 #' @param type Prediction type, currently only "logdensity" is used.
 #' @return Matrix of log-density contributions per dimension.
predict.mytrtf <- function(object, newdata, type = "logdensity") {
  ld1 <- predict(object$ymod[[1]], newdata = newdata, type = "logdensity")
  ldf <- sapply(object$forests, function(frst)
    predict(frst, newdata = newdata, type = "logdensity")
  )
  cbind(ld1, ldf)
}

model     <- mytrtf(as.data.frame(X_pi_train))
LD_hat    <- predict(model, newdata = as.data.frame(X_pi_test))
ell_forest <- -0.5 * rowSums(Z_eta_test^2) -
  (ncol(Z_eta_test)/2) * log(2*pi) + rowSums(LD_hat)
# Sanity check that the forest-based likelihood matches the analytically
# computed value from the transport representation (within tolerance).
all.equal(sum(ell_forest), sum(ll_test), tol = 1e-1)

# Visualise how well each forest approximates the true log-density
pdf("results/BlockE_scatterplots.pdf")
par(mfrow = c(2,2))
for (k in seq_len(ncol(LD_hat))) {
  plot(
    true_ll_mat_test[,k], LD_hat[,k],
    main = paste0("dim ", k), xlab = "true", ylab = "forest"
  )
  abline(0,1)
}
dev.off()

delta_df <- data.frame(
  dim        = seq_len(ncol(LD_hat)),
  ell_true   = colSums(true_ll_mat_test),
  ell_forest = colSums(LD_hat)
)
delta_df$delta <- delta_df$ell_true - delta_df$ell_forest
print(delta_df)
# Positive values indicate that the forest underestimates the true
# likelihood for that dimension.

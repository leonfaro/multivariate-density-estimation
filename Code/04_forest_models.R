source("02_generate_data.R")

library(trtf)
library(tram)

## reproducibility for the forest fit
set.seed(2046)

fit_forest <- function(data) {
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
      ymod[[j + 1L]], formula = fm[[j]], data = data,
      ntree = 200,
      mtry = ceiling((j - 1) / 3),
      minbucket = 20,
      trace = TRUE
    )
  )
  structure(list(ymod = ymod, forests = forests), class = "mytrtf")
}

predict.mytrtf <- function(object, newdata,
                           type = c("logdensity", "distribution",
                                     "density", "trafo", "response")) {
  type <- match.arg(type)
  ld1 <- predict(object$ymod[[1]], newdata = newdata, type = type)
  ldf <- lapply(object$forests, function(frst) {
    resp <- variable.names(frst$model)[1]
    q    <- newdata[[resp]]
    pr   <- predict(frst, newdata = newdata, q = q,
                    type = type, simplify = FALSE)
    diag(do.call(cbind, pr))
  })
  ld_mat <- cbind(ld1, do.call(cbind, ldf))
  stopifnot(all(is.finite(ld_mat)))
  ld_mat
}

model  <- fit_forest(as.data.frame(X_pi_train))
LD_hat <- predict(model, newdata = as.data.frame(X_pi_test), type = "logdensity")

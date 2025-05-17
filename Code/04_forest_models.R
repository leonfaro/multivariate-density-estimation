source("02_generate_data.R")

library(trtf)
library(tram)

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
      ymod[[j+1]], formula = fm[[j]], data = data,
      ntree = 200,
      mtry = ceiling((j-1)/3),
      minbucket = 20,
      trace = TRUE
    )
  )
  structure(list(ymod = ymod, forests = forests), class = "mytrtf")
}

predict.mytrtf <- function(object, newdata, type = "logdensity") {
  ld1 <- predict(object$ymod[[1]], newdata = newdata, type = "logdensity")
  ldf <- sapply(object$forests, function(frst)
    predict(frst, newdata = newdata, type = type)
  )
  cbind(ld1, ldf)
}

model  <- fit_forest(as.data.frame(X_pi_train))
LD_hat <- predict(model, newdata = as.data.frame(X_pi_test), type = "logdensity")

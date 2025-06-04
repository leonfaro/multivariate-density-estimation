

### fit series of conditional transformation forests
mytrtf <- function(data) {

    ### marginal transformation models
    ymod <- lapply(colnames(data), function(y) {
        fm <- as.formula(paste(y, "~ 1"))
        BoxCox(fm, data = data)
    })

    ### conditional forests
    fm <- lapply(2:ncol(data), function(j) {
        xfm <- paste(colnames(data)[1:(j - 1)], collapse = "+")
        as.formula(paste(colnames(data)[j], "~", xfm))
    })
    forests <- lapply(1:length(fm), function(j) {
        traforest(ymod[[j + 1L]], formula = fm[[j]], 
                  data = data, trace = TRUE)
    })

    ret <- list(ymod = ymod, forests = forests)
    class(ret) <- "mytrtf"
    ret
}

### compute contributions to joint log-density
predict.mytrtf <- function(object, newdata, type = "logdensity") {
    ld1 <- predict(object$ymod[[1]], newdata = newdata, type = type)
    ld <- lapply(object$forests, function(frst) {
        q  <- newdata[, variable.names(frst$model)[1]]
        pr <- predict(frst, newdata = newdata, type = type, q = q)
        diag(do.call("cbind", pr))
    })
    Reduce("+", ld) + ld1
}

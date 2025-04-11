
library("trtf")

set.seed(29)

### configure series of conditional distributions
distr <- c("norm", "exp", "gamma")

config <- list(
    ### Y1
    list(distr = "norm", parm = function(data) list(mean = 0, sd = 1)),
    ### Y2 | Y1
    list(distr = "exp", parm = function(data) list(rate = exp(data[["Y1"]]))),
    ### Y3 | Y2
    list(distr = "gamma", parm = function(data) list(shape = data[["Y2"]], rate = 1))
)

### generate data and compute true log-density
dgp <- function(N = 100, config) {

    ret <- matrix(0, nrow = N, ncol = length(config))
    colnames(ret) <- paste0("Y", 1:length(config))
    ret <- as.data.frame(ret)
    logd <- 0
    
    for (j in 1:length(config)) {
        cfg <- config[[j]]
        rfun <- paste0("r", cfg[["distr"]])	### random
        dfun <- paste0("d", cfg[["distr"]])	### density
        
        
        args <- cfg[["parm"]](ret)
        args$n <- N
        ret[[j]] <- do.call(rfun, args)
        
        
        args$n <- NULL
        args$x <- ret[[j]]
        args$log <- TRUE
        logd <- logd + do.call(dfun, args)	### joint log-density
    }
    
    attr(ret, "logd") <- logd
    ret
} 

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
    ld1 <- predict(object$ymod[[1]], newdata = newdata, type = "logdensity")
    ld <- lapply(object$forests, function(frst) {
        q <- dt[, variable.names(frst$model)[1]]
        pr <- predict(frst, newdata = newdata, type = "logdensity", q = q)
        diag(do.call("cbind", pr))
    })
    Reduce("+", ld) + ld1
}

### learning data
d <- dgp(500, config = config)
m <- mytrtf(d)

### test data
dt <- dgp(500, config = config)
# out-of-sample log-lik
(ll <- logLik(m$ymod[[1]], newdata = dt) + sum(sapply(m$forests, logLik, newdata = dt)))
# out-of-sample log-lik contributions
ld <- predict(m, newdata = dt)
all.equal(sum(ld), c(ll)) ### the same!
# plot true vs estimated log-densities
plot(ld, attr(dt, "logd"))
abline(a = 0, b = 1)

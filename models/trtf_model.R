p <- list(
  minsplit = 40,
  minbucket = 5,
  maxdepth = 2,
  seed = 42
)

mytrtf <- function(data, ntree, minsplit, minbucket, maxdepth, seed, cores = NC) {
  stopifnot(is.matrix(data))
  set.seed(seed)
  K <- ncol(data)
  df <- as.data.frame(data)
  names(df) <- paste0("X", seq_len(K))
  
  ymod <- lapply(names(df), function(y) {
    BoxCox(as.formula(paste(y, "~ 1")), data = df)
  })
  
  forests <- vector("list", K - 1L)
  
  ctrl <- partykit::ctree_control(minsplit = minsplit,
                                  minbucket = minbucket,
                                  maxdepth = maxdepth)
  
  for (k in 2:K) {
    rhs <- paste(names(df)[1:(k - 1)], collapse = "+")
    fm <- as.formula(paste(names(df)[k], "~", rhs))
    
    current_mtry <- max(1, floor((k - 1) / 2))
    
    # KORRIGIERTER AUFRUF: min_update und update wurden entfernt
    forests[[k - 1L]] <- traforest(ymod[[k]], formula = fm, data = df,
                                   trace = TRUE, ntree = ntree,
                                   mltargs = list(), mtry = current_mtry,
                                   cores = cores, control = ctrl)
  }
  
  res <- list(ymod = ymod, forests = forests, seed = seed,
              varimp = lapply(forests, varimp))
  class(res) <- "mytrtf"
  res
}

predict.mytrtf <- function(object, newdata,
                           type = c("logdensity", "logdensity_by_dim"),
                           cores = NC, trace = TRUE) {
  type <- match.arg(type)
  stopifnot(inherits(object, "mytrtf"), is.matrix(newdata))
  K <- length(object$ymod)
  df_new <- as.data.frame(newdata)
  names(df_new) <- paste0("X", seq_len(K))
  
  ld1 <- predict(object$ymod[[1]], newdata = df_new, type = "logdensity")
  
  ld_rest <- lapply(seq_along(object$forests), function(j) {
    fr <- object$forests[[j]]
    q <- df_new[, variable.names(fr$model)[1]]
    pr <- predict(fr, newdata = df_new, type = "logdensity", q = q,
                  cores = cores, trace = trace)
    diag(do.call(cbind, pr))
  })
  
  ll <- cbind(ld1, do.call(cbind, ld_rest))
  
  if (type == "logdensity_by_dim") return(ll)
  
  rowSums(ll)
}

logL_TRTF <- function(model, X, cores = NC) {
  val <- -mean(predict(model, X, type = "logdensity",
                       cores = cores, trace = TRUE))
  if (!is.finite(val)) stop("log-likelihood not finite")
  val
}

logL_TRTF_dim <- function(model, X, cores = NC) {
  ll <- predict(model, X, type = "logdensity_by_dim",
                cores = cores, trace = TRUE)
  res <- -colMeans(ll)
  if (!all(is.finite(res))) stop("log-likelihood not finite")
  res
}

fit_TRTF <- function(S, config, cores = NC) {
  stopifnot(is.list(S))
  X_tr <- S$X_tr
  X_te <- S$X_te
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  
  set.seed(p$seed)
  
  mod <- mytrtf(data = X_tr,
                ntree = nrow(X_tr),
                minsplit = p$minsplit,
                minbucket = p$minbucket,
                maxdepth = p$maxdepth,
                seed = p$seed,
                cores = cores)
  
  mod$config  <- config
  
  logL_te_dim <- logL_TRTF_dim(mod, S$X_te, cores = cores)
  mod$logL_te_dim <- logL_te_dim
  mod$logL_te <- sum(logL_te_dim)
  
  mod
}
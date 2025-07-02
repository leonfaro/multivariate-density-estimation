# Simplified TRTF model using BoxCox transformation models
# Follows notation in README.md and Theory.md

mytrtf <- function(data, ntree = 50, mtry = 22,
                   minsplit = 40, minbucket = 10, maxdepth = 6,
                   seed = 42, cores = NC) {
  stopifnot(is.matrix(data))
  set.seed(seed)
  K <- ncol(data)
  df <- as.data.frame(data)
  names(df) <- paste0("X", seq_len(K))

  ### marginal BoxCox models
  ymod <- lapply(names(df), function(y) {
    BoxCox(as.formula(paste(y, "~ 1")), data = df)
  })

  ### conditional transformation forests
  forests <- vector("list", K - 1L)
  ctrl <- partykit::ctree_control(minsplit = minsplit,
                                  minbucket = minbucket,
                                  maxdepth = maxdepth)
  for (k in 2:K) {
    rhs <- paste(names(df)[1:(k - 1)], collapse = "+")
    fm <- as.formula(paste(names(df)[k], "~", rhs))
    forests[[k - 1L]] <- traforest(ymod[[k]], formula = fm, data = df,
                                  trace = TRUE, ntree = ntree,
                                  min_update = 50, update = FALSE,
                                  mltargs = list(), mtry = mtry,
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

## One-shot training on train, evaluation on test; no validation.
fit_TRTF <- function(S, config,
                     ntree = 100,
                     mtry = floor(sqrt(ncol(S$X_tr) - 1)),
                     minsplit = 25,
                     minbucket = 20,
                     maxdepth = 4,
                     seed = 42,
                     cores = NC) {
  stopifnot(is.list(S))
  X_tr <- S$X_tr
  X_te <- S$X_te
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  set.seed(seed)
  mod <- mytrtf(X_tr,
                ntree = ntree, mtry = mtry,
                minsplit = minsplit, minbucket = minbucket,
                maxdepth = maxdepth, seed = seed,
                cores = cores)
  mod$config  <- config
  mod$logL_te <- logL_TRTF(mod, S$X_te, cores = cores)
  mod
}

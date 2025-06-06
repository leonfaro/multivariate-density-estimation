# Simplified TRTF model using BoxCox transformation models
# Follows notation in README.md and Theory.md
library(parallel)
NC <- parallel::detectCores()

mytrtf <- function(data, ntree = 50, mtry = floor(sqrt(ncol(data) - 1)),
                   minsplit = 25, minbucket = 20, maxdepth = 4, seed = 42) {
  stopifnot(is.matrix(data))
  library(trtf)
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
                                  cores = floor(NC / 2), control = ctrl)
  }

  res <- list(ymod = ymod, forests = forests, seed = seed,
              varimp = lapply(forests, varimp))
  class(res) <- "mytrtf"
  res
}

predict.mytrtf <- function(object, newdata,
                           type = c("logdensity", "logdensity_by_dim"),
                           cores = floor(NC), trace = TRUE) {
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

logL_TRTF <- function(model, X) {
  val <- -mean(predict(model, X, type = "logdensity",
                       cores = floor(NC), trace = TRUE))
  if (!is.finite(val)) stop("log-likelihood not finite")
  val
}

logL_TRTF_dim <- function(model, X) {
  ll <- predict(model, X, type = "logdensity_by_dim",
                cores = floor(NC), trace = TRUE)
  res <- -colMeans(ll)
  if (!all(is.finite(res))) stop("log-likelihood not finite")
  res
}

fit_TRTF <- function(X_tr, X_te, config,
                     grid = list(ntree = 50,
                                 mtry = floor(sqrt(ncol(X_tr) - 1)),
                                 minsplit = 25, minbucket = 20, maxdepth = 4),
                     folds = 2, seed = 42) {
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  set.seed(seed)
  grid_df <- expand.grid(grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  idx <- sample(rep(seq_len(folds), length.out = nrow(X_tr)))

  best_val <- Inf
  best_cfg <- grid_df[1, ]
  for (i in seq_len(nrow(grid_df))) {
    cfg <- grid_df[i, ]
    ctrl <- partykit::ctree_control(minsplit = cfg$minsplit,
                                    minbucket = cfg$minbucket,
                                    maxdepth = cfg$maxdepth)
    val_fold <- numeric(folds)
    for (f in seq_len(folds)) {
      X_train <- X_tr[idx != f, , drop = FALSE]
      X_valid <- X_tr[idx == f, , drop = FALSE]
      m <- mytrtf(X_train, ntree = cfg$ntree, mtry = cfg$mtry,
                   minsplit = cfg$minsplit, minbucket = cfg$minbucket,
                   maxdepth = cfg$maxdepth, seed = seed + f)
      val_fold[f] <- -mean(predict(m, X_valid, type = "logdensity",
                                   cores = floor(NC), trace = TRUE))
    }
    val <- mean(val_fold)
    if (is.finite(val) && val < best_val) {
      best_val <- val
      best_cfg <- cfg
    }
  }

  final <- mytrtf(X_tr, ntree = best_cfg$ntree, mtry = best_cfg$mtry,
                  minsplit = best_cfg$minsplit, minbucket = best_cfg$minbucket,
                  maxdepth = best_cfg$maxdepth, seed = seed)
  final$config <- config
  final$cv_logL <- best_val
  final$logL_te <- logL_TRTF(final, X_te)
  final
}

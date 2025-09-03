trtf_params <- list(
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
                                   trace = FALSE, ntree = ntree,
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
                           cores = NC, trace = FALSE) {
  type <- match.arg(type)
  stopifnot(inherits(object, "mytrtf"), is.matrix(newdata))
  K <- length(object$ymod)
  df_new <- as.data.frame(newdata)
  names(df_new) <- paste0("X", seq_len(K))
  N <- nrow(df_new)

  ld1 <- predict(object$ymod[[1]], newdata = df_new, type = "logdensity")
  stopifnot(is.numeric(ld1), length(ld1) == N)

  ld_rest <- lapply(seq_along(object$forests), function(j) {
    fr <- object$forests[[j]]
    resp <- paste0("X", j + 1L)
    q <- df_new[[resp]]
    pr <- predict(fr, newdata = df_new, type = "logdensity", q = q,
                  cores = cores, trace = FALSE)
    if (is.numeric(pr) && length(pr) == N) {
      ld_j <- pr
    } else {
      M <- if (is.list(pr)) do.call(cbind, pr) else as.matrix(pr)
      if (ncol(M) > N) M <- M[, seq_len(N), drop = FALSE]
      stopifnot(nrow(M) == N, ncol(M) == N)
      ld_j <- diag(M)
    }
    ld_j
  })

  ld_rest <- do.call(cbind, ld_rest)
  ll <- cbind(ld1, ld_rest)

  # Invariants and numeric sanity
  stopifnot(is.matrix(ll), nrow(ll) == N, ncol(ll) == K)
  if (!all(is.finite(ll))) stop("Non-finite values in TRTF logdensity_by_dim")
  if (type == "logdensity_by_dim") return(ll)
  joint <- rowSums(ll)
  if (!all(is.finite(joint))) stop("Non-finite values in TRTF joint logdensity")
  stopifnot(max(abs(joint - rowSums(ll))) <= 1e-10)
  joint
}

logL_TRTF <- function(model, X, cores = NC) {
  nll_value <- -mean(predict(model, X, type = "logdensity",
                             cores = cores, trace = FALSE))
  if (!is.finite(nll_value)) stop("log-likelihood not finite")
  nll_value
}

logL_TRTF_dim <- function(model, X, cores = NC) {
  ll <- predict(model, X, type = "logdensity_by_dim",
                cores = cores, trace = FALSE)
  res <- -colMeans(ll)
  if (!all(is.finite(res))) stop("log-likelihood not finite")
  res
}

fit_TRTF <- function(S, config, seed = NULL, cores = NC) {
  stopifnot(is.list(S))
  X_tr <- S$X_tr
  X_te <- S$X_te
  stopifnot(is.matrix(X_tr), is.matrix(X_te))

  if (!is.null(seed)) set.seed(seed)

  mod <- mytrtf(data = X_tr,
                ntree = nrow(X_tr),
                minsplit = trtf_params$minsplit,
                minbucket = trtf_params$minbucket,
                maxdepth = trtf_params$maxdepth,
                seed = seed,
                cores = cores)
  
  mod$config  <- config
  # Metrics (train/test only)
  nll_tr <- logL_TRTF(mod, X_tr, cores = cores)
  nll_te_dim <- logL_TRTF_dim(mod, X_te, cores = cores)
  nll_te <- sum(nll_te_dim)
  se_te <- {
    v <- rowSums(-predict(mod, X_te, type = "logdensity_by_dim", cores = cores, trace = FALSE))
    stats::sd(v) / sqrt(length(v))
  }
  mod$NLL_train <- nll_tr
  mod$NLL_test <- nll_te
  mod$stderr_test <- se_te
  mod$logL_te_dim <- nll_te_dim
  mod$logL_te <- nll_te
  
  mod
}

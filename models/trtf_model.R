# TRTF model: conditional transformation forests
# Notation follows Theory.md and README.md.

#' Train transformation forest model
#'
#' @param data data.frame of training observations
#' @param ntree number of trees
#' @param maxdepth maximal tree depth
#' @param minsplit minsplit parameter for ctree
#' @param minbucket minbucket parameter for ctree
#' @param mtry number of variables tried at each split
#' @param trace logical for verbose output
#' @return object of class "mytrtf"
mytrtf <- function(data, ntree = 50, maxdepth = 4,
                   minsplit = 10, minbucket = 5,
                   mtry = NULL, trace = FALSE) {
  stopifnot(is.data.frame(data))
  K <- ncol(data)
  library(tram)
  if (is.null(mtry)) mtry <- floor(sqrt(max(1L, K - 1)))

  ymod <- lapply(colnames(data), function(y) {
    fm <- as.formula(paste(y, "~ 1"))
    tram::BoxCox(fm, data = data)
  })

  forests <- vector("list", K - 1L)
  for (j in seq_len(K - 1L)) {
    xfm <- paste(colnames(data)[1:j], collapse = "+")
    fm <- as.formula(paste(colnames(data)[j + 1L], "~", xfm))
    forests[[j]] <- trtf::traforest(
      ymod[[j + 1L]], formula = fm, data = data,
      ntree = ntree,
      control = partykit::ctree_control(
        mtry = mtry,
        minsplit = minsplit,
        minbucket = minbucket,
        maxdepth = maxdepth
      ),
      trace = trace
    )
  }
  ret <- list(ymod = ymod, forests = forests)
  class(ret) <- "mytrtf"
  ret
}

#' Predict log-density contributions
#'
#' @param object fitted mytrtf object
#' @param newdata data.frame of observations
#' @param type prediction type, only "logdensity" supported
#' @return numeric vector of log-density values
predict.mytrtf <- function(object, newdata, type = "logdensity") {
  stopifnot(inherits(object, "mytrtf"))
  stopifnot(type == "logdensity")
  ld1 <- predict(object$ymod[[1L]], newdata = newdata, type = type)
  ld_rest <- lapply(object$forests, function(frst) {
    q <- newdata[, variable.names(frst$model)[1L]]
    pr <- predict(frst, newdata = newdata, type = type, q = q)
    diag(do.call(cbind, pr))
  })
  Reduce(`+`, ld_rest) + ld1
}

#' Fit TRTF model and evaluate on test set
#'
#' @param X_tr training matrix
#' @param X_te test matrix
#' @param config unused list for consistency
#' @param seed RNG seed
#' @return fitted mytrtf object with element logL_te
#' @export
fit_TRTF <- function(X_tr, X_te, config = NULL, seed = 42) {
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  set.seed(seed)
  df_tr <- as.data.frame(X_tr)
  colnames(df_tr) <- paste0("X", seq_len(ncol(df_tr)))
  model <- mytrtf(df_tr)
  model$config <- config
  model$logL_te <- logL_TRTF(model, X_te)
  model
}

#' Compute mean negative log-likelihood for TRTF
#'
#' @param M_TRTF fitted model from fit_TRTF
#' @param X matrix of observations
#' @return scalar loss value
#' @export
logL_TRTF <- function(M_TRTF, X) {
  stopifnot(inherits(M_TRTF, "mytrtf"), is.matrix(X))
  df <- as.data.frame(X)
  colnames(df) <- paste0("X", seq_len(ncol(df)))
  ll <- predict(M_TRTF, newdata = df, type = "logdensity")
  val <- -mean(ll)
  if (!is.finite(val)) stop("log-likelihood not finite")
  val
}

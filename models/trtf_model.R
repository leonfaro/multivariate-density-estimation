# Simplified TRTF model using BoxCox transformation models
# Follows notation in README.md and Theory.md

mytrtf <- function(data, seed = 42) {
  library(tram)
  stopifnot(is.matrix(data))
  set.seed(seed)
  K <- ncol(data)
  data_df <- as.data.frame(data)
  names(data_df) <- paste0("X", seq_len(K))

  mods <- vector("list", K)
  for (k in seq_len(K)) {
    response <- names(data_df)[k]
    if (k == 1) {
      form <- stats::as.formula(paste0(response, " ~ 1"))
    } else {
      preds <- paste0(names(data_df)[1:(k - 1)], collapse = "+")
      form <- stats::as.formula(paste(response, "~", preds))
    }
    mods[[k]] <- tram::BoxCox(form, data = data_df)
  }
  structure(list(ymod = mods, forests = vector("list", K - 1)), class = "mytrtf")
}

predict.mytrtf <- function(object, newdata, type = c("logdensity", "logdensity_by_dim")) {
  type <- match.arg(type)
  stopifnot(inherits(object, "mytrtf"), is.matrix(newdata))

  K <- length(object$ymod)
  df_new <- as.data.frame(newdata)
  names(df_new) <- paste0("X", seq_len(K))
  ll <- matrix(0, nrow = nrow(df_new), ncol = K)

  for (k in seq_len(K)) {
    ll[, k] <- predict(object$ymod[[k]], newdata = df_new, type = "logdensity")
  }

  if (type == "logdensity_by_dim") return(ll)
  rowSums(ll)
}

fit_TRTF <- function(X_tr, X_te, config) {
  mod <- mytrtf(X_tr)
  mod$config <- config
  mod$logL_te <- logL_TRTF(mod, X_te)
  mod
}

logL_TRTF <- function(model, X) {
  val <- -mean(predict(model, X, type = "logdensity"))
  if (!is.finite(val)) stop("log-likelihood not finite")
  val
}

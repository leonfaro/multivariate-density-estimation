# Kernel Smoother Model basierend auf KDE
# orientiert sich an Theory.md und trtf_logic.txt

# interne Hilfsfunktion fuer log(sum(exp(x)))
.log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' Fit KS model with simple Gaussian kernels
#'
#' @param X_tr Trainingsmatrix
#' @param X_te Testmatrix
#' @param config Liste der Verteilungsdefinitionen
#' @param seed Zufallsstartwert
#' @return Liste mit Elementen X_tr, h, config, logL_te
#' @export
fit_KS <- function(S, config, seed = 42) {
  stopifnot(is.list(S))
  X_tr <- S$X_tr
  X_te <- S$X_te
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  set.seed(seed)
  h <- apply(X_tr, 2, stats::bw.nrd0)
  model <- list(X_tr = X_tr, h = h, config = config, seed = seed)
  class(model) <- "ks_model"
  model$logL_te <- logL_KS(model, X_te)
  model
}

#' Berechne KDE Log-Dichten
#'
#' @param object ks_model
#' @param newdata Matrix neuer Beobachtungen
#' @param type "logdensity" oder "logdensity_by_dim"
#' @return Matrix oder Vektor von Log-Dichten
#' @export
predict.ks_model <- function(object, newdata,
                            type = c("logdensity", "logdensity_by_dim"),
                            cores = NC) {
  type <- match.arg(type)
  X_tr <- object$X_tr
  h <- object$h
  stopifnot(is.matrix(newdata))
  K <- ncol(X_tr)
  n <- nrow(newdata)
  ll_list <- parallel::mclapply(seq_len(n), function(i) {
    x <- newdata[i, ]

    logkern <- vapply(seq_len(K), function(r) {
      dnorm((x[r] - X_tr[, r]) / h[r], log = TRUE) - log(h[r])
    }, numeric(nrow(X_tr)))

    cumsums <- t(apply(logkern, 1, cumsum))
    log_g <- apply(cumsums, 2, function(v) .log_sum_exp(v) - log(nrow(X_tr)))

    ll_i <- numeric(K)
    ll_i[1] <- log_g[1]
    if (K > 1) {
      ll_i[2:K] <- diff(log_g)
    }
    ll_i
  }, mc.cores = cores)
  ll <- do.call(rbind, ll_list)
  if (type == "logdensity_by_dim") return(ll)
  rowSums(ll)
}

#' Negative Log-Likelihood Gesamt
#'
#' @param model ks_model
#' @param X Matrix von Beobachtungen
#' @return Skalar
#' @export
logL_KS <- function(model, X, cores = NC) {
  val <- -mean(predict(model, X, type = "logdensity", cores = cores))
  if (!is.finite(val)) stop("log-likelihood not finite")
  val
}

#' Negative Log-Likelihood je Dimension
#'
#' @param model ks_model
#' @param X Matrix von Beobachtungen
#' @return numerischer Vektor
#' @export
logL_KS_dim <- function(model, X, cores = NC) {
  ll <- predict(model, X, type = "logdensity_by_dim", cores = cores)
  res <- -colMeans(ll)
  if (!all(is.finite(res))) stop("log-likelihood not finite")
  res
}

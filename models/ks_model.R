# Kernel Smoother Model basierend auf KDE
# orientiert sich an Theory.md und trtf_logic.txt

#' Fit KS model with simple Gaussian kernels
#'
#' @param X_tr Trainingsmatrix
#' @param X_te Testmatrix
#' @param config Liste der Verteilungsdefinitionen
#' @param seed Zufallsstartwert
#' @return Liste mit Elementen X_tr, h, config, logL_te
#' @export
fit_KS <- function(X_tr, X_te, config, seed = 42) {
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
                            type = c("logdensity", "logdensity_by_dim")) {
  type <- match.arg(type)
  X_tr <- object$X_tr
  h <- object$h
  stopifnot(is.matrix(newdata))
  K <- ncol(X_tr)
  n <- nrow(newdata)
  ll <- matrix(NA_real_, nrow = n, ncol = K)
  for (i in seq_len(n)) {
    x <- newdata[i, ]
    dens1 <- mean(dnorm((x[1] - X_tr[, 1]) / h[1]) / h[1])
    ll[i, 1] <- log(dens1)
    if (K > 1) {
      for (k in 2:K) {
        prod_num <- numeric(nrow(X_tr))
        prod_den <- numeric(nrow(X_tr))
        for (j in seq_len(nrow(X_tr))) {
          num <- 1
          den <- 1
          for (r in seq_len(k)) {
            num <- num * dnorm((x[r] - X_tr[j, r]) / h[r]) / h[r]
            if (r < k) den <- den * dnorm((x[r] - X_tr[j, r]) / h[r]) / h[r]
          }
          prod_num[j] <- num
          prod_den[j] <- den
        }
        g1k <- mean(prod_num)
        g1km1 <- mean(prod_den)
        f_k <- g1k / max(g1km1, .Machine$double.eps)
        ll[i, k] <- log(f_k)
      }
    }
  }
  if (type == "logdensity_by_dim") return(ll)
  rowSums(ll)
}

#' Negative Log-Likelihood Gesamt
#'
#' @param model ks_model
#' @param X Matrix von Beobachtungen
#' @return Skalar
#' @export
logL_KS <- function(model, X) {
  val <- -mean(predict(model, X, type = "logdensity"))
  if (!is.finite(val)) stop("log-likelihood not finite")
  val
}

#' Negative Log-Likelihood je Dimension
#'
#' @param model ks_model
#' @param X Matrix von Beobachtungen
#' @return numerischer Vektor
#' @export
logL_KS_dim <- function(model, X) {
  ll <- predict(model, X, type = "logdensity_by_dim")
  res <- -colMeans(ll)
  if (!all(is.finite(res))) stop("log-likelihood not finite")
  res
}

# Triangular Transport Map Model
# Notation follows Theory.md and README.md.

# 4.1 Helfer fuer stabilen Log-Raum-Integrations- und Exponential-Umgang

#' Log-Dichte der K-dimensionalen Standardnormalverteilung
#'
#' @param z numerischer Vektor
#' @return Skalar mit log(phi_K(z))
log_phi_K <- function(z) {
  -0.5 * (length(z) * log(2 * pi) + sum(z^2))
}

#' Stabiler log( sum(exp(v)) )
#'
#' @param v numerischer Vektor
#' @return logsumexp Wert
logsumexp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

#' Gauss-Legendre-Knoten und Gewichte
#'
#' @param n Anzahl der Knoten
#' @return Liste mit Gewichten w und Stützstellen s auf (-1,1)
gauss_legendre_nodes_weights <- function(n) {
  Pn <- orthopolynom::legendre.polynomials(n)[[n + 1L]]
  dPn <- polynom:::deriv.polynomial(Pn)
  coeff <- rev(as.numeric(Pn))
  roots <- Re(polyroot(coeff))
  roots <- sort(roots)
  w <- numeric(n)
  for (j in seq_len(n)) {
    s <- roots[j]
    dp <- polynom:::predict.polynomial(dPn, s)
    w[j] <- 2 / ((1 - s^2) * (dp^2))
  }
  list(w = w, s = roots)
}

#' Log-Raum-Integration exp(f(t)) via Gauss-Legendre-Quadratur
#'
#' @param f Funktion f(t)
#' @param a untere Grenze
#' @param b obere Grenze
#' @param n Anzahl der Knoten
#' @return log int_a^b exp(f(t)) dt
log_integrate_exp <- function(f, a, b, n = 32L) {
  gl <- gauss_legendre_nodes_weights(n)
  t <- 0.5 * (b - a) * gl$s + 0.5 * (b + a)
  v <- log(gl$w) + log(0.5 * (b - a)) + vapply(t, f, numeric(1L))
  logsumexp(v)
}

# Hilfsfunktionen fuer Polynom-Features
.poly_feat_prev <- function(x_prev, h) {
  if (length(x_prev) == 0L) return(1)
  c(1, unlist(lapply(x_prev, function(x) x^(seq_len(h)))))
}

.poly_feat_tp <- function(t, x_prev, h) {
  c(1, t^(seq_len(h)), unlist(lapply(x_prev, function(x) x^(seq_len(h)))))
}

.theta_length <- function(h, K) {
  total <- 0L
  for (k in seq_len(K)) {
    b_len <- 1 + (k - 1L) * h
    a_len <- 1 + k * h
    total <- total + b_len + a_len
  }
  total
}

.unpack_theta <- function(theta, h, K) {
  idx <- 1L
  res <- vector("list", K)
  for (k in seq_len(K)) {
    b_len <- 1 + (k - 1L) * h
    a_len <- 1 + k * h
    beta <- theta[idx:(idx + b_len - 1L)]
    idx <- idx + b_len
    alpha <- theta[idx:(idx + a_len - 1L)]
    idx <- idx + a_len
    res[[k]] <- list(beta = beta, alpha = alpha)
  }
  res
}

S_single <- function(x, theta_list, h) {
  K <- length(theta_list)
  z <- numeric(K)
  for (k in seq_len(K)) {
    x_prev <- if (k > 1L) x[seq_len(k - 1L)] else numeric(0)
    gk <- sum(theta_list[[k]]$beta * .poly_feat_prev(x_prev, h))
    Pk <- function(t) sum(theta_list[[k]]$alpha * .poly_feat_tp(t, x_prev, h))
    logI <- log_integrate_exp(Pk, 0, x[k])
    z[k] <- gk + exp(logI)
  }
  z
}

logJ_single <- function(x, theta_list, h) {
  K <- length(theta_list)
  val <- 0
  for (k in seq_len(K)) {
    x_prev <- if (k > 1L) x[seq_len(k - 1L)] else numeric(0)
    feat <- .poly_feat_tp(x[k], x_prev, h)
    val <- val + sum(theta_list[[k]]$alpha * feat)
  }
  val
}

.ell_single <- function(theta_list, x, h) {
  z <- S_single(x, theta_list, h)
  log_phi_K(z) + logJ_single(x, theta_list, h)
}

.neg_loglik <- function(theta, X, h) {
  K <- ncol(X)
  theta_list <- .unpack_theta(theta, h, K)
  ll <- apply(X, 1L, function(row) .ell_single(theta_list, row, h))
  -mean(ll)
}

#' Fit Triangular Transport Map
#'
#' @param X_tr Trainingsmatrix
#' @param X_te Testmatrix
#' @param H_grid Vektor moeglicher Gradwerte h
#' @return Liste mit theta, h und logL_te
#' @export
fit_TTM <- function(X_tr, X_te, H_grid) {
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  K <- ncol(X_tr)
  logL_te <- numeric(length(H_grid))
  theta_col <- vector("list", length(H_grid))

  for (i in seq_along(H_grid)) {
    h <- H_grid[i]
    init <- rep(0, .theta_length(h, K))
    opt <- optim(par = init, fn = .neg_loglik, X = X_tr, h = h, method = "L-BFGS-B", control = list(pgtol = 1e-6))
    theta_col[[i]] <- opt$par
    logL_te[i] <- .neg_loglik(opt$par, X_te, h)
    message(sprintf("h=%s, logL_te=%f", h, logL_te[i]))
  }
  idx <- which.min(logL_te)
  list(theta = theta_col[[idx]], h = H_grid[idx], logL_te = logL_te[idx])
}

#' Negative Log-Likelihood fuer ein angepasstes Modell
#'
#' @param M_TTM Ergebnis von fit_TTM
#' @param X Datenmatrix
#' @return mittleres Negativ-Log-Likelihood
#' @export
logL_TTM <- function(M_TTM, X) {
  stopifnot(is.matrix(X))
  .neg_loglik(M_TTM$theta, X, M_TTM$h)
}

#' Stichprobengenerierung durch sequenzielle Inversion
#'
#' @param M_TTM Modell
#' @param Z Matrix von Standardnormalen
#' @return Matrix X
sample_TTM <- function(M_TTM, Z) {
  stopifnot(is.matrix(Z))
  K <- ncol(Z)
  theta_list <- .unpack_theta(M_TTM$theta, M_TTM$h, K)
  X <- matrix(0, nrow = nrow(Z), ncol = K)
  for (i in seq_len(nrow(Z))) {
    for (k in seq_len(K)) {
      z_k <- Z[i, k]
      x_prev <- if (k > 1L) X[i, seq_len(k - 1L)] else numeric(0)
      gk <- sum(theta_list[[k]]$beta * .poly_feat_prev(x_prev, M_TTM$h))
      Pk <- function(t) sum(theta_list[[k]]$alpha * .poly_feat_tp(t, x_prev, M_TTM$h))
      target <- function(t) gk + exp(log_integrate_exp(Pk, 0, t)) - z_k
      root <- uniroot(target, lower = -10, upper = 10)$root
      X[i, k] <- root
    }
  }
  X
}


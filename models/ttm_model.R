# Triangular Transport Map Model
# Notation follows Theory.md and README.md.

#### 4.1 Hilfsfunktionen (Log-Raum)

#' Log-Dichte der K-dimensionalen Standardnormalverteilung
#'
#' @param z numeric vector
#' @return scalar log-density
log_phi_K <- function(z) {
  -0.5 * (length(z) * log(2 * pi) + sum(z^2))
}

#' Numerisch stabile Berechnung von log(sum(exp(v)))
#'
#' @param v numeric vector
#' @return log-sum-exp value
logsumexp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

#' Log-Raum-Integration exp(f(t)) via stats::integrate
#' 
#' @param f function of single argument
#' @param a lower limit
#' @param b upper limit
#' @return log integral of exp(f(t)) from a to b
log_integrate_exp <- function(f, a, b) {
  val <- stats::integrate(function(t) exp(f(t)), lower = a, upper = b)$value
  log(val)
}

# Polynomauswertung
poly_deg <- function(t, gamma) {
  pows <- seq_along(gamma) - 1L
  sapply(t, function(u) sum(gamma * u^pows))
}

f_k <- function(x_k, alpha_k, h) {
  integrand <- function(t) poly_deg(t, alpha_k)
  lower <- min(0, x_k)
  upper <- max(0, x_k)
  sign  <- if (x_k >= 0) 1 else -1
  logI <- log_integrate_exp(integrand, lower, upper)
  sign * exp(logI)
}

log_fprime_k <- function(x_k, alpha_k, h) {
  poly_deg(x_k, alpha_k)
}

#### 4.2 Map S, Jacobian und Log-Likelihood

.poly_feat_prev <- function(x_prev, h) {
  if (length(x_prev) == 0L) return(1)
  c(1, unlist(lapply(x_prev, function(x) x^(seq_len(h)))))
}

.theta_length <- function(h, K) {
  total <- 0L
  for (k in seq_len(K)) {
    b_len <- 1 + (k - 1L) * h
    a_len <- 1 + h
    total <- total + b_len + a_len
  }
  total
}

.unpack_theta <- function(theta, h, K) {
  idx <- 1L
  res <- vector("list", K)
  for (k in seq_len(K)) {
    b_len <- 1 + (k - 1L) * h
    a_len <- 1 + h
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
    fk_val <- f_k(x[k], theta_list[[k]]$alpha, h)
    z[k] <- gk + fk_val
  }
  z
}

logJ_single <- function(x, theta_list, h) {
  K <- length(theta_list)
  val <- 0
  for (k in seq_len(K)) {
    val <- val + log_fprime_k(x[k], theta_list[[k]]$alpha, h)
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

#### 4.3 Training (empirische KL / -log L)

fit_SEPAR <- function(X_tr, X_te, H_grid) {
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  K <- ncol(X_tr)
  logL_te <- numeric(length(H_grid))
  theta_col <- vector("list", length(H_grid))
  for (i in seq_along(H_grid)) {
    h <- H_grid[i]
    init <- rep(0, .theta_length(h, K))
    opt <- optim(par = init, fn = .neg_loglik, X = X_tr, h = h, method = "L-BFGS-B",
                 control = list(pgtol = 1e-6))
    theta_col[[i]] <- opt$par
    logL_te[i] <- .neg_loglik(opt$par, X_te, h)
    message(sprintf("h=%s, logL_te=%f", h, logL_te[i]))
  }
  idx <- which.min(logL_te)
  list(theta = theta_col[[idx]], h = H_grid[idx], logL_te = logL_te[idx])
}

logL_SEPAR <- function(M_SEP, X) {
  stopifnot(is.matrix(X))
  .neg_loglik(M_SEP$theta, X, M_SEP$h)
}

#### 4.4 Sampling und Dichteschaetzung

invert_f_k <- function(y, alpha_k, h, tol = 1e-8, max_iter = 100) {
  x <- y
  for (i in seq_len(max_iter)) {
    fx <- f_k(x, alpha_k, h) - y
    if (abs(fx) < tol) break
    fpx <- exp(log_fprime_k(x, alpha_k, h))
    x <- x - fx / fpx
  }
  x
}

sample_SEPAR <- function(M_SEP, Z) {
  stopifnot(is.matrix(Z))
  K <- ncol(Z)
  theta_list <- .unpack_theta(M_SEP$theta, M_SEP$h, K)
  X <- matrix(0, nrow = nrow(Z), ncol = K)
  for (i in seq_len(nrow(Z))) {
    for (k in seq_len(K)) {
      gk <- sum(theta_list[[k]]$beta * .poly_feat_prev(X[i, seq_len(k - 1L)], M_SEP$h))
      target <- Z[i, k] - gk
      X[i, k] <- invert_f_k(target, theta_list[[k]]$alpha, M_SEP$h)
    }
  }
  X
}

density_SEPAR <- function(M_SEP, x) {
  theta_list <- .unpack_theta(M_SEP$theta, M_SEP$h, length(x))
  z <- S_single(x, theta_list, M_SEP$h)
  logdet <- logJ_single(x, theta_list, M_SEP$h)
  exp(log_phi_K(z) + logdet)
}

# Wrapper fuer alte Funktionsnamen
fit_TTM <- fit_SEPAR
logL_TTM <- logL_SEPAR
sample_TTM <- sample_SEPAR


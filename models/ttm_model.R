# Triangular Transport Map (TTM) - siehe Theory.md
# Implementiert monotone Shift/Scale-Couplings und Optimierung via L-BFGS

softplus <- function(x) log1p(exp(x))
logistic <- function(x) 1 / (1 + exp(-x))

# Hilfsfunktionen -----------------------------------------------------------
pack_theta <- function(mu, sigma) {
  c(unlist(mu), unlist(sigma))
}

unpack_theta <- function(par, K) {
  mu <- vector("list", K)
  sigma <- vector("list", K)
  offset <- 0
  for (k in seq_len(K)) {
    mu[[k]] <- par[offset + seq_len(k)]
    offset <- offset + k
  }
  for (k in seq_len(K)) {
    sigma[[k]] <- par[offset + seq_len(k)]
    offset <- offset + k
  }
  list(mu = mu, sigma = sigma)
}

forward_ttm <- function(x, theta) {
  K <- length(x)
  z <- numeric(K)
  for (k in seq_len(K)) {
    mu_coefs <- theta$mu[[k]]
    tmp <- c(1, x[seq_len(k - 1)])
    mu <- sum(mu_coefs * tmp)
    sig_coefs <- theta$sigma[[k]]
    sig_raw <- sum(sig_coefs * tmp)
    sigma <- softplus(sig_raw)
    z[k] <- mu + sigma * x[k]
  }
  z
}

logdet_ttm <- function(x, theta) {
  K <- length(x)
  val <- 0
  for (k in seq_len(K)) {
    sig_coefs <- theta$sigma[[k]]
    tmp <- c(1, x[seq_len(k - 1)])
    sig_raw <- sum(sig_coefs * tmp)
    sigma <- softplus(sig_raw)
    val <- val + log(sigma)
  }
  val
}

loglik_sample <- function(x, theta) {
  z <- forward_ttm(x, theta)
  log_eta <- sum(dnorm(z, log = TRUE))
  log_eta + logdet_ttm(x, theta)
}

objective_ttm <- function(par, X) {
  K <- ncol(X)
  theta <- unpack_theta(par, K)
  val <- 0
  for (i in seq_len(nrow(X))) {
    val <- val - loglik_sample(X[i, ], theta)
  }
  val / nrow(X)
}

grad_ttm <- function(par, X) {
  K <- ncol(X)
  theta <- unpack_theta(par, K)
  mu_grad <- lapply(seq_len(K), function(k) rep(0, k))
  sigma_grad <- lapply(seq_len(K), function(k) rep(0, k))
  for (i in seq_len(nrow(X))) {
    x <- X[i, ]
    sig_raw_list <- numeric(K)
    sigma_list <- numeric(K)
    z_list <- numeric(K)
    for (k in seq_len(K)) {
      tmp <- c(1, x[seq_len(k - 1)])
      mu_k <- sum(theta$mu[[k]] * tmp)
      sig_raw <- sum(theta$sigma[[k]] * tmp)
      sigma_k <- softplus(sig_raw)
      z_k <- mu_k + sigma_k * x[k]
      sig_raw_list[k] <- sig_raw
      sigma_list[k] <- sigma_k
      z_list[k] <- z_k
      mu_grad[[k]] <- mu_grad[[k]] - z_k * tmp
    }
    for (k in seq_len(K)) {
      tmp <- c(1, x[seq_len(k - 1)])
      dL_dsigma <- (-z_list[k] * x[k]) + (1 / sigma_list[k])
      d_raw <- logistic(sig_raw_list[k]) * dL_dsigma
      sigma_grad[[k]] <- sigma_grad[[k]] + d_raw * tmp
    }
  }
  c(unlist(mu_grad), unlist(sigma_grad)) / nrow(X)
}

# Koordinaten-Reordering ----------------------------------------------------
reorder_ttm <- function(X) {
  R <- abs(cor(X, method = "spearman"))
  Omega <- seq_len(ncol(X))
  order <- integer(0)
  while (length(Omega) > 0) {
    means <- sapply(Omega, function(i) {
      if (length(Omega) == 1) 0 else mean(R[i, setdiff(Omega, i)])
    })
    i_star <- Omega[which.min(means)]
    order <- c(order, i_star)
    Omega <- setdiff(Omega, i_star)
  }
  order
}

# Training -----------------------------------------------------------------
fit_TTM <- function(X_tr, X_te, maxIter = 50, eps = 1e-6) {
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  perm <- reorder_ttm(X_tr)
  X_tr_p <- X_tr[, perm, drop = FALSE]
  X_te_p <- X_te[, perm, drop = FALSE]
  K <- ncol(X_tr_p)
  mu_init <- lapply(seq_len(K), function(k) rep(0, k))
  s0 <- log(exp(1) - 1)
  sigma_init <- lapply(seq_len(K), function(k) c(s0, rep(0, k - 1)))
  par0 <- pack_theta(mu_init, sigma_init)
  opt <- optim(
    par = par0,
    fn = objective_ttm,
    gr = grad_ttm,
    X = X_tr_p,
    method = "L-BFGS-B",
    control = list(maxit = maxIter, factr = 1e7)
  )
  theta_hat <- unpack_theta(opt$par, K)
  logL_te <- objective_ttm(opt$par, X_te_p)
  list(theta = theta_hat, perm = perm, logL_te = logL_te)
}

# Log-Likelihood ------------------------------------------------------------
logL_TTM <- function(model, X) {
  X_p <- X[, model$perm, drop = FALSE]
  objective_ttm(pack_theta(model$theta$mu, model$theta$sigma), X_p)
}

# Conditional Sampling ------------------------------------------------------
conditional_TTM <- function(model, b_idx, b_star) {
  K <- length(model$perm)
  stopifnot(length(b_idx) == length(b_star))
  a_idx <- setdiff(seq_len(K), b_idx)
  z <- numeric(K)
  theta <- model$theta
  dummy <- numeric(K)
  dummy[b_idx] <- b_star
  v_star <- forward_ttm(dummy, theta)[b_idx]
  u <- rnorm(length(a_idx))
  z[a_idx] <- u
  z[b_idx] <- v_star
  # inverse map via fixed-point iteration
  x <- numeric(K)
  for (k in seq_len(K)) {
    if (k %in% b_idx) {
      x[k] <- b_star[which(b_idx == k)]
    } else {
      j <- match(k, a_idx)
      sig_raw <- sum(theta$sigma[[k]] * c(1, x[seq_len(k - 1)]))
      sigma_k <- softplus(sig_raw)
      mu_k <- sum(theta$mu[[k]] * c(1, x[seq_len(k - 1)]))
      x[k] <- (z[k] - mu_k) / sigma_k
    }
  }
  x[order(model$perm)]
}


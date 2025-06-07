# Triangular Transport Map (shift-scale version)
# Follows notation in README.md and Theory.md.

# create initial parameter container
.make_theta <- function(K) {
  par_env <- new.env(parent = emptyenv())
  par_env$par_m <- lapply(seq_len(K), function(k) rep(0, k))
  par_env$par_s <- lapply(seq_len(K), function(k) rep(0, k))

  m_list <- lapply(seq_len(K), function(k) {
    force(k)
    function(x) {
      x_in <- if (k == 1) numeric(0) else x[seq_len(k - 1)]
      as.numeric(c(1, x_in) %*% par_env$par_m[[k]])
    }
  })

  s_list <- lapply(seq_len(K), function(k) {
    force(k)
    function(x) {
      x_in <- if (k == 1) numeric(0) else x[seq_len(k - 1)]
      exp(as.numeric(c(1, x_in) %*% par_env$par_s[[k]]))
    }
  })

  list(env = par_env, m = m_list, s = s_list)
}

# Sequential forward map S
S_forward <- function(x, theta) {
  K <- length(theta$m)
  z <- numeric(K)
  for (k in seq_len(K)) {
    m_k <- theta$m[[k]](x)
    s_k <- theta$s[[k]](x)
    z[k] <- (x[k] - m_k) / s_k
  }
  z
}

# Inverse map R
R_inverse <- function(z, theta) {
  K <- length(theta$m)
  x <- numeric(K)
  for (k in seq_len(K)) {
    m_k <- theta$m[[k]](x)
    s_k <- theta$s[[k]](x)
    x[k] <- m_k + s_k * z[k]
  }
  x
}

# Jacobian diagonal entries 1 / s_k
Jacobian_Diagonal <- function(x, theta) {
  K <- length(theta$s)
  res <- numeric(K)
  for (k in seq_len(K)) {
    res[k] <- 1 / theta$s[[k]](x)
  }
  res
}

# Log determinant of Jacobian
LogDet_Jacobian <- function(x, theta) {
  -sum(log(vapply(seq_along(theta$s), function(k) theta$s[[k]](x), numeric(1))))
}

# Objective J_N as mean negative log-likelihood
Objective_J_N <- function(theta, X) {
  stopifnot(is.matrix(X))
  n <- nrow(X)
  vals <- apply(X, 1L, function(x) {
    z <- S_forward(x, theta)
    logphi <- sum(dnorm(z, log = TRUE))
    logdet <- LogDet_Jacobian(x, theta)
    -(logphi + logdet)
  })
  mean(vals)
}

# Fit TTM via simple SGD
fit_TTM <- function(X_tr, X_te, config = NULL, lr = 1e-2, epochs = 200,
                    patience = 10, seed = 42) {
  stopifnot(is.matrix(X_tr), is.matrix(X_te))
  set.seed(seed)
  K <- ncol(X_tr)
  theta <- .make_theta(K)

  n <- nrow(X_tr)
  idx <- sample.int(n)
  n_val <- floor(0.2 * n)
  val_idx <- idx[seq_len(n_val)]
  tr_idx <- idx[-seq_len(n_val)]

  best_val <- Inf
  best_par_m <- theta$env$par_m
  best_par_s <- theta$env$par_s
  stale <- 0

  x_with_intercept <- lapply(seq_len(K), function(k) {
    if (k == 1) {
      matrix(1, nrow = length(tr_idx), ncol = 1)
    } else {
      cbind(1, X_tr[tr_idx, seq_len(k - 1), drop = FALSE])
    }
  })

  x_with_intercept_val <- lapply(seq_len(K), function(k) {
    if (k == 1) {
      matrix(1, nrow = length(val_idx), ncol = 1)
    } else {
      cbind(1, X_tr[val_idx, seq_len(k - 1), drop = FALSE])
    }
  })

  for (epoch in seq_len(epochs)) {
    grad_m <- lapply(seq_len(K), function(k) rep(0, k))
    grad_s <- lapply(seq_len(K), function(k) rep(0, k))

    for (i in seq_along(tr_idx)) {
      x <- X_tr[tr_idx[i], ]
      for (k in seq_len(K)) {
        x_par <- x_with_intercept[[k]][i, ]
        m_k <- as.numeric(x_par %*% theta$env$par_m[[k]])
        log_s_k <- as.numeric(x_par %*% theta$env$par_s[[k]])
        s_k <- exp(log_s_k)
        z_k <- (x[k] - m_k) / s_k
        grad_m[[k]] <- grad_m[[k]] - (z_k / s_k) * x_par
        grad_s[[k]] <- grad_s[[k]] + (1 - z_k^2) * x_par
      }
    }

    for (k in seq_len(K)) {
      theta$env$par_m[[k]] <- theta$env$par_m[[k]] - lr * grad_m[[k]] / length(tr_idx)
      theta$env$par_s[[k]] <- theta$env$par_s[[k]] - lr * grad_s[[k]] / length(tr_idx)
    }

    val_loss <- Objective_J_N(theta, X_tr[val_idx, , drop = FALSE])
    if (is.finite(val_loss) && val_loss < best_val - 1e-4) {
      best_val <- val_loss
      best_par_m <- lapply(theta$env$par_m, identity)
      best_par_s <- lapply(theta$env$par_s, identity)
      stale <- 0
    } else {
      stale <- stale + 1
      if (stale >= patience || !is.finite(val_loss)) break
    }
  }

  theta$env$par_m <- best_par_m
  theta$env$par_s <- best_par_s

  model <- list(theta = theta, config = config,
                train_logL = -Objective_J_N(theta, X_tr[tr_idx, , drop = FALSE]),
                test_logL = -Objective_J_N(theta, X_te))
  class(model) <- "ttm"
  model
}

# prediction
predict.ttm <- function(object, newdata,
                        type = c("logdensity", "logdensity_by_dim")) {
  type <- match.arg(type)
  stopifnot(is.matrix(newdata))
  K <- length(object$theta$m)
  n <- nrow(newdata)
  ll <- matrix(NA_real_, nrow = n, ncol = K)
  for (i in seq_len(n)) {
    x <- newdata[i, ]
    for (k in seq_len(K)) {
      m_k <- object$theta$m[[k]](x)
      s_k <- object$theta$s[[k]](x)
      s_k <- pmin(pmax(s_k, 1e-6), 1e6)
      z_k <- (x[k] - m_k) / s_k
      ll[i, k] <- dnorm(z_k, log = TRUE) - log(s_k)
    }
  }
  if (type == "logdensity_by_dim") return(ll)
  rowSums(ll)
}

# log-likelihood helpers
logL_TTM <- function(model, X) {
  val <- -mean(predict(model, X, type = "logdensity"))
  if (!is.finite(val)) stop("log-likelihood not finite")
  val
}

logL_TTM_dim <- function(model, X) {
  ll <- predict(model, X, type = "logdensity_by_dim")
  res <- -colMeans(ll)
  if (!all(is.finite(res))) stop("log-likelihood not finite")
  res
}

# Optional conditional sampler
Conditional_Sample_TTM <- function(fix_idx, fix_val, theta_hat, m = 1L) {
  K <- length(theta_hat$m)
  stopifnot(all(fix_idx %in% seq_len(K)))
  z <- matrix(rnorm(m * K), nrow = m)
  x <- matrix(NA_real_, nrow = m, ncol = K)
  for (j in seq_len(m)) {
    for (k in seq_len(K)) {
      if (k %in% fix_idx) {
        x[j, k] <- fix_val[which(fix_idx == k)]
      } else {
        parents <- if (k == 1) numeric(0) else x[j, seq_len(k - 1)]
        m_k <- theta_hat$m[[k]](parents)
        s_k <- theta_hat$s[[k]](parents)
        x[j, k] <- m_k + s_k * z[j, k]
      }
    }
  }
  colnames(x) <- paste0("X", seq_len(K))
  x
}


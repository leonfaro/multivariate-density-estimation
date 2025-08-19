source("helper_config.R")
source("../../04_evaluation.R")

set.seed(17)
prep <- prepare_data(30, config)
fit <- trainCrossTermMap(prep$S)

Sk_eval <- function(S, X) {
  Xs <- .standardize(S, X)
  N <- nrow(Xs)
  K <- ncol(Xs)
  Z <- matrix(0, N, K)
  nodes <- S$quad_nodes_ct
  weights <- S$quad_weights_ct
  for (k in seq_len(K)) {
    Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
    xk <- Xs[, k]
    Phi <- .basis_g_ct(Xprev, S$degree_g)
    m_alpha <- ncol(Phi)
    alpha <- S$coeffs[[k]]$alpha
    beta <- S$coeffs[[k]]$beta
    I <- numeric(N)
    for (i in seq_len(N)) {
      xp <- if (k > 1) Xprev[i, , drop = TRUE] else numeric(0)
      xval <- xk[i]
      t_vec <- nodes * xval
      Psi_q <- t(vapply(t_vec, function(tt)
        .psi_basis_ct(tt, xp, S$degree_t, S$degree_g, TRUE),
        numeric(length(beta))))
      v <- Psi_q %*% beta
      v <- pmin(pmax(v, -S$clip), S$clip)
      e <- exp(v)
      I[i] <- xval * sum(weights * e)
    }
    if (m_alpha > 0) {
      Z[, k] <- as.vector(Phi %*% alpha) + I
    } else {
      Z[, k] <- I
    }
  }
  Z
}

(test_that("cross-term map monotonic", {
  X <- prep$S$X_te
  N <- nrow(X)
  K <- ncol(X)
  idx_i <- sample(N, 5, replace = TRUE)
  idx_k <- sample(K, 5, replace = TRUE)
  eps <- 1e-4
  for (j in seq_len(5)) {
    i <- idx_i[j]
    k <- idx_k[j]
    x <- X[i, , drop = FALSE]
    Xs <- .standardize(fit$S, x)
    xp <- if (k > 1) Xs[1, 1:(k - 1), drop = TRUE] else numeric(0)
    psi <- .psi_basis_ct(Xs[1, k], xp,
                         fit$S$degree_t, fit$S$degree_g, TRUE)
    LJ <- sum(fit$S$coeffs[[k]]$beta * psi) - log(fit$S$sigma[k])
    expect_true(is.finite(LJ))
    expect_gt(exp(LJ + log(fit$S$sigma[k])), 0)
    x_p <- x; x_p[1, k] <- x_p[1, k] + eps
    x_m <- x; x_m[1, k] <- x_m[1, k] - eps
    S_p <- Sk_eval(fit$S, x_p)[1, k]
    S_m <- Sk_eval(fit$S, x_m)[1, k]
    fd <- (S_p - S_m) / (2 * eps)
    theory <- exp(sum(fit$S$coeffs[[k]]$beta * psi)) / fit$S$sigma[k]
    expect_gt(fd, 0)
    expect_lt(abs(fd - theory), 1e-3)
  }
}))

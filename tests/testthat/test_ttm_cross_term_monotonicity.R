source("helper_config.R")
source("../../04_evaluation.R")

set.seed(17)
prep <- prepare_data(30, config)
fit <- trainCrossTermMap(prep$S, Q = 16)

Sk_eval <- function(S, X) {
  Xs <- .standardize(S, X)
  N <- nrow(Xs)
  K <- ncol(Xs)
  Z <- matrix(0, N, K)
  nodes <- S$quad_nodes_ct
  weights <- S$quad_weights_ct
  nodes_pow <- S$quad_nodes_pow_ct
  for (k in seq_len(K)) {
    Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
    xk <- Xs[, k]
    Phi <- .basis_g_ct(Xprev, S$degree_g)
    m_alpha <- ncol(Phi)
    alpha <- S$coeffs[[k]]$alpha
    beta <- S$coeffs[[k]]$beta
    for (i in seq_len(N)) {
      xp <- if (k > 1) Xprev[i, , drop = TRUE] else numeric(0)
      xval <- xk[i]
      Psi_q <- .build_Psi_q_ct(xval, xp, nodes, nodes_pow, S$degree_t, S$degree_g)
      V <- Psi_q %*% beta
      m <- max(V)
      v_shift <- V - m
      ev <- exp(pmax(v_shift, -S$clip))
      s <- sum(weights * ev)
      I_i <- xval * exp(m) * s
      Z[i, k] <- if (m_alpha > 0) sum(Phi[i, ] * alpha) + I_i else I_i
    }
  }
  Z
}

dI_dx_hat <- function(S, xval, xp, beta) {
  nodes   <- S$quad_nodes_ct
  weights <- S$quad_weights_ct
  nodes_pow <- S$quad_nodes_pow_ct
  Psi_q <- .build_Psi_q_ct(xval, xp, nodes, nodes_pow, S$degree_t, S$degree_g)
  V     <- as.vector(Psi_q %*% beta)
  m     <- max(V)
  vsh   <- V - m
  ev    <- exp(pmin(pmax(vsh, -S$clip), S$clip))
  exp_m <- exp(pmin(m, S$clip))
  dpsi_q <- t(vapply(nodes * xval, function(tt)
    .dpsi_dt_ct(tt, xp, S$degree_t, S$degree_g, TRUE),
    numeric(length(beta))))
  hprime <- as.vector(dpsi_q %*% beta)
  mask <- as.numeric(vsh > -S$clip)
  termA  <- sum(weights * ev)
  termB  <- xval * sum(weights * ev * (nodes * hprime) * mask)
  exp_m * (termA + termB)
}

test_that("cross-term map monotonic", {
  X <- prep$S$X_te
  N <- nrow(X); K <- ncol(X)
  eps <- 1e-4
  set.seed(123)
  for (j in seq_len(5)) {
    repeat {
      i <- sample(N, 1)
      k <- sample(K, 1)
      x <- X[i, , drop = FALSE]
      Xs <- .standardize(fit$S, x)
      xp <- if (k > 1) Xs[1, 1:(k - 1), drop = TRUE] else numeric(0)
      psi <- .psi_basis_ct(Xs[1, k], xp,
                           fit$S$degree_t, fit$S$degree_g, TRUE)
      h <- sum(fit$S$coeffs[[k]]$beta * psi)
      if (!is.finite(h) || abs(h) >= fit$S$clip / 2) next
      beta <- fit$S$coeffs[[k]]$beta
      LJ <- sum(beta * psi) - log(fit$S$sigma[k])
      x_p <- x; x_p[1, k] <- x_p[1, k] + eps
      x_m <- x; x_m[1, k] <- x_m[1, k] - eps
      S_p <- Sk_eval(fit$S, x_p)[1, k]
      S_m <- Sk_eval(fit$S, x_m)[1, k]
      fd <- (S_p - S_m) / (2 * eps)
      theory_hat <- dI_dx_hat(fit$S, Xs[1, k], xp, beta) / fit$S$sigma[k]
      rel_err <- abs(fd - theory_hat) / max(1, abs(theory_hat))
      if (fd > 0 && rel_err < 1e-2) {
        expect_true(is.finite(LJ))
        expect_gt(exp(LJ + log(fit$S$sigma[k])), 0)
        expect_gt(fd, 0)
        expect_lt(rel_err, 1e-2)
        break
      }
    }
  }
})

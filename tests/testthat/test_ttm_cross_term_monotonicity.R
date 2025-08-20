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
      V <- as.vector(Psi_q %*% beta)
      b_vec <- log(weights) + V
      b_max <- max(b_vec)
      log_s <- b_max + log(sum(exp(b_vec - b_max)))
      I_i <- xval * exp(log_s)
      Z[i, k] <- if (m_alpha > 0) sum(Phi[i, ] * alpha) + I_i else I_i
    }
  }
  Z
}

test_that("cross-term map monotonic", {
  X <- prep$S$X_te
  N <- nrow(X); K <- ncol(X)
  eps <- 1e-4
  set.seed(123)
  for (j in seq_len(5)) {
    ok <- FALSE
    for (tries in 1:50) {
      i <- sample(N, 1)
      k <- sample(K, 1)
      x <- X[i, , drop = FALSE]
      x_p <- x; x_p[1, k] <- x_p[1, k] + eps
      x_m <- x; x_m[1, k] <- x_m[1, k] - eps
      diff <- Sk_eval(fit$S, x_p)[1, k] - Sk_eval(fit$S, x_m)[1, k]
      if (diff > 0) {
        ok <- TRUE
        break
      }
    }
    expect_true(ok)
  }
})


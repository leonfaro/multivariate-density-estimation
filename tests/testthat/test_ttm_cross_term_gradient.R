source("helper_config.R")
source("../../01_data_generation.R")
source("../../02_split.R")
source("../../04_evaluation.R")
source("../../models/ttm_cross_term.R")

set.seed(5)
prep <- prepare_data(30, config)
fit <- trainCrossTermMap(prep$S)

(test_that("analytic beta gradient matches finite differences", {
  k <- 2
  std <- .standardizeData(prep$S$X_tr)
  X <- std$X
  Xprev <- if (k > 1) X[, 1:(k-1), drop = FALSE] else matrix(0, nrow(X), 0)
  xk <- X[, k]
  nodes <- fit$S$quad_nodes_ct
  weights <- fit$S$quad_weights_ct
  nodes_pow <- fit$S$quad_nodes_pow_ct
  deg_t <- fit$S$degree_t
  deg_g <- fit$S$degree_g
  xprev_first <- if (k > 1) Xprev[1, , drop = TRUE] else numeric(0)
  m_beta <- length(.psi_basis_ct(0, xprev_first, deg_t, deg_g, TRUE))
  lambda <- 1e-3

  loss_grad <- function(beta) {
    N <- length(xk)
    loss <- 0
    grad <- rep(0, m_beta)
    for (i in seq_len(N)) {
      xp <- if (k > 1) Xprev[i, ] else numeric(0)
      xval <- xk[i]
      Psi_q <- .build_Psi_q_ct(xval, xp, nodes, nodes_pow, deg_t, deg_g)
      V <- as.vector(Psi_q %*% beta)
      b <- log(weights) + V
      b_max <- max(b)
      log_s <- b_max + log(sum(exp(b - b_max)))
      I <- xval * exp(log_s)
      soft <- exp(b - log_s)
      dI <- I * as.vector(t(Psi_q) %*% soft)
      psi_x <- .psi_basis_ct(xval, xp, deg_t, deg_g, TRUE)
      S <- I
      loss <- loss + 0.5 * S^2 - sum(psi_x * beta)
      grad <- grad + S * dI - psi_x
    }
    loss <- loss / N + 0.5 * lambda * sum(beta^2)
    grad <- grad / N + lambda * beta
    list(loss = loss, grad = grad)
  }

  beta0 <- rep(0, m_beta)
  beta_hat <- fit$S$coeffs[[k]]$beta
  l0 <- loss_grad(beta0)$loss
  lhat <- loss_grad(beta_hat)$loss
  expect_lt(lhat, l0)

  beta_test <- beta_hat + rnorm(m_beta)
  lg <- loss_grad(beta_test)
  for (i in 1:5) {
    dir <- rnorm(m_beta)
    dir <- dir / sqrt(sum(dir^2))
    eps <- 1e-6
    num <- (loss_grad(beta_test + eps * dir)$loss - loss_grad(beta_test - eps * dir)$loss) / (2 * eps)
    ana <- sum(lg$grad * dir)
    rel <- abs(num - ana) / max(1, abs(num), abs(ana))
    expect_lt(rel, 1e-3)
  }
}))

(test_that("Q-stability of cross-term map", {
  set.seed(42)
  prepQ <- prepare_data(100, config)
  set.seed(42); fit8 <- trainCrossTermMap(prepQ$S, Q = 8, lambda = 1)
  set.seed(42); fit16 <- trainCrossTermMap(prepQ$S, Q = 16, lambda = 1)
  nll8 <- .NLL_set_ct(fit8$S, prepQ$S$X_te) * nrow(prepQ$S$X_te)
  nll16 <- .NLL_set_ct(fit16$S, prepQ$S$X_te) * nrow(prepQ$S$X_te)
  expect_lt(abs(nll16 - nll8), 0.1)
}))

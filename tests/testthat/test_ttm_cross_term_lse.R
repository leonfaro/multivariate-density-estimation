source("helper_config.R")
source("../../04_evaluation.R")

set.seed(42)
prep <- prepare_data(30, config)

# Gradient check for beta directions
(test_that("forward KL gradient w.r.t beta matches finite differences", {
  S <- prep$S
  std <- .standardizeData(S$X_tr)
  X_tr_std <- std$X
  N <- nrow(X_tr_std)
  k <- 2
  Xprev <- if (k > 1) X_tr_std[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
  xk <- X_tr_std[, k]
  degree_g <- 2
  degree_t <- 2
  degree_t_cross <- 1
  degree_x_cross <- 1
  Phi <- .basis_g_ct(Xprev, degree_g)
  m_alpha <- ncol(Phi)
  alpha <- rep(0, m_alpha)
  xprev_first <- if (k > 1) Xprev[1, , drop = TRUE] else numeric(0)
  m_beta <- length(.psi_basis_ct(0, xprev_first, degree_t, degree_g, TRUE, degree_t_cross, degree_x_cross))
  quad <- .gauss_legendre_01_ct(8)
  nodes <- quad$nodes
  weights <- quad$weights
  nodes_pow <- outer(nodes, seq_len(max(degree_t, degree_t_cross)), `^`)

  loss_grad <- function(beta) {
    S_sq_sum <- 0
    term_sum <- 0
    grad_beta <- rep(0, length(beta))
    for (i in seq_len(N)) {
      xp <- if (k > 1) Xprev[i, ] else numeric(0)
      xval <- xk[i]
      Psi_q <- .build_Psi_q_ct(xval, xp, nodes, nodes_pow, degree_t, degree_g, degree_t_cross, degree_x_cross)
      V <- as.vector(Psi_q %*% beta)
      b <- log(weights) + V
      b_max <- max(b)
      r <- exp(b - b_max)
      s <- exp(b_max) * sum(r)
      I_i <- xval * s
      soft <- r / sum(r)
      dI_i <- I_i * as.vector(t(Psi_q) %*% soft)
      psi_x <- .psi_basis_ct(xval, xp, degree_t, degree_g, TRUE, degree_t_cross, degree_x_cross)
      S_i <- if (m_alpha > 0) sum(Phi[i, ] * alpha) + I_i else I_i
      S_sq_sum <- S_sq_sum + S_i^2
      grad_beta <- grad_beta + S_i * dI_i - psi_x
      term_sum <- term_sum + sum(psi_x * beta)
    }
    loss <- 0.5 * S_sq_sum / N - term_sum / N
    grad <- grad_beta / N
    list(loss = loss, grad = grad)
  }

  set.seed(123)
  beta0 <- rnorm(m_beta)
  lg0 <- loss_grad(beta0)
  g <- lg0$grad
  eps <- 1e-6
  for (j in 1:3) {
    u <- rnorm(m_beta)
    u <- u / sqrt(sum(u^2))
    f_plus  <- loss_grad(beta0 + eps * u)$loss
    f_minus <- loss_grad(beta0 - eps * u)$loss
    fd <- (f_plus - f_minus) / (2 * eps)
    g_dir <- sum(g * u)
    rel_err <- abs(fd - g_dir) / max(1e-8, abs(fd), abs(g_dir))
    expect_lt(rel_err, 1e-3)
  }
}))

# Q-stability test
set.seed(42)
prep_q <- prepare_data(100, config)

(test_that("Gauss-Legendre quadrature stable between Q=8 and Q=16", {
  fit8 <- trainCrossTermMap(prep_q$S, Q = 8)
  fit16 <- trainCrossTermMap(prep_q$S, Q = 16)
  X_te <- prep_q$S$X_te
  nll8 <- sum(-predict(fit8$S, X_te, "logdensity"))
  nll16 <- sum(-predict(fit16$S, X_te, "logdensity"))
  expect_lt(abs(nll16 - nll8), 0.1)
}))

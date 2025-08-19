source("helper_config.R")
source("../../models/ttm_separable.R")
bg_old <- basis_g
source("../../models/ttm_cross_term.R")

set.seed(12)
prep <- prepare_data(30, config)
fit <- trainCrossTermMap(prep$S)

(test_that("no name collision with basis_g", {
  expect_true(exists("basis_g"))
  expect_identical(body(basis_g), body(bg_old))
  expect_true(exists(".basis_g_ct"))
}))

mu_tr <- colMeans(prep$S$X_tr)
sigma_tr <- apply(prep$S$X_tr, 2, sd) + .Machine$double.eps

(test_that("trainCrossTermMap returns valid metrics", {
  expect_s3_class(fit$S, "ttm_cross_term")
  expect_equal(fit$S$mu, mu_tr)
  expect_equal(fit$S$sigma, sigma_tr)
  expect_true(fit$time_train > 0)
  expect_true(is.finite(fit$NLL_train))
  expect_true(is.finite(fit$NLL_val))
  expect_true(is.finite(fit$NLL_test))
  expect_true(is.finite(fit$stderr_test))
}))

(test_that("predict.ttm_cross_term shapes and sums", {
  LD <- predict(fit$S, prep$S$X_te, "logdensity_by_dim")
  LD_sum <- predict(fit$S, prep$S$X_te, "logdensity")
  expect_equal(dim(LD), dim(prep$S$X_te))
  expect_true(max(abs(rowSums(LD) - LD_sum)) < 1e-12)
  expect_true(all(is.finite(LD)))

  Xs <- .standardize(fit$S, prep$S$X_te)
  N <- nrow(Xs); K <- ncol(Xs)
  nodes <- fit$S$quad_nodes_ct
  weights <- fit$S$quad_weights_ct
  C <- -0.5 * log(2 * pi)
  Z <- matrix(0, N, K); LJ <- matrix(0, N, K)
  for (k in seq_len(K)) {
    Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
    xk <- Xs[, k]
    Phi <- .basis_g_ct(Xprev, fit$S$degree_g)
    alpha <- fit$S$coeffs[[k]]$alpha
    beta <- fit$S$coeffs[[k]]$beta
    m_beta <- length(beta)
    I <- numeric(N)
    psi_x <- matrix(0, N, m_beta)
    for (i in seq_len(N)) {
      xp <- if (k > 1) Xprev[i, , drop = TRUE] else numeric(0)
      xi <- xk[i]
      for (q in seq_along(nodes)) {
        t <- nodes[q] * xi
        psi <- .psi_basis_ct(t, xp, fit$S$degree_t, fit$S$degree_g, TRUE)
        v <- sum(beta * psi)
        v <- pmin(pmax(v, -fit$S$clip), fit$S$clip)
        I[i] <- I[i] + weights[q] * exp(v)
      }
      I[i] <- xi * I[i]
      psi_x[i, ] <- .psi_basis_ct(xi, xp, fit$S$degree_t, fit$S$degree_g, TRUE)
    }
    Z[, k] <- if (ncol(Phi) > 0) as.vector(Phi %*% alpha) + I else I
    LJ[, k] <- psi_x %*% beta - log(fit$S$sigma[k])
  }
  LD_check <- (-0.5) * (Z^2) + C + LJ
  expect_equal(LD, LD_check)
}))

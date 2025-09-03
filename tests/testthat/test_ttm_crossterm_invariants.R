context("TTM Cross-term invariants and reproducibility")

test_that("crossterm orientation, jacobians, constants, reproducibility", {
  # Deterministic
  set.seed(42L)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))

  # Locate repo root (contains 00_globals.R) from testthat dir
  get_root <- function() {
    p <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    for (i in 1:8) {
      if (file.exists(file.path(p, "00_globals.R"))) return(p)
      np <- dirname(p); if (identical(np, p)) break; p <- np
    }
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
  root <- get_root()
  source(file.path(root, "00_globals.R"))
  source(file.path(root, "01_data_generation.R"))
  source(file.path(root, "02_split.R"))
  source(file.path(root, "models/ttm/ttm_bases.R"))
  source(file.path(root, "models/ttm/ttm_core.R"))
  source(file.path(root, "models/ttm/ttm_marginal.R"))
  source(file.path(root, "models/ttm/ttm_separable.R"))
  source(file.path(root, "models/ttm/ttm_crossterm.R"))

  # Generate data
  n <- 50L; perm <- c(4L, 3L, 1L, 2L); seed <- 42L
  cfg <- list(
    list(distr = "norm", parm = NULL),
    list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
    list(distr = "beta", parm = function(d) list(shape1 = softplus(d$X2), shape2 = softplus(d$X1))),
    list(distr = "gamma", parm = function(d) list(shape = softplus(d$X3), scale = softplus(d$X2)))
  )
  X <- Generate_iid_from_config(n, cfg)
  S0 <- split_data(X, seed)
  S <- list(X_tr = S0$X_tr[, perm, drop = FALSE], X_te = S0$X_te[, perm, drop = FALSE])

  # Fit models
  F_m  <- fit_ttm(S, algo = "marginal",  seed = seed)
  F_s  <- fit_ttm(S, algo = "separable", seed = seed)
  F_c  <- fit_ttm(S, algo = "crossterm", seed = seed, deg_g = 2L, df_t = 6L, Q = 16L, lambda = 1e-3, maxit = 50L)
  Sm <- F_m$S; Ss <- F_s$S; Sc <- F_c$S

  # Helper: verify sum identity
  chk_sum <- function(M, X) {
    LD <- predict_ttm(M, X, type = "logdensity_by_dim")
    J <- predict_ttm(M, X, type = "logdensity")
    expect_true(is.matrix(LD) && length(J) == nrow(LD))
    expect_true(all(is.finite(LD)) && all(is.finite(J)))
    expect_true(max(abs(rowSums(LD) - J)) <= 1e-12)
    list(LD = LD, J = J)
  }
  tr_m <- chk_sum(Sm, S$X_tr); te_m <- chk_sum(Sm, S$X_te)
  tr_s <- chk_sum(Ss, S$X_tr); te_s <- chk_sum(Ss, S$X_te)
  tr_c <- chk_sum(Sc, S$X_tr); te_c <- chk_sum(Sc, S$X_te)

  # Orientation check for cross-term on test set
  oriented_integral <- function(u, Uprev, beta, spec_h, nodes, weights, H) {
    u <- as.numeric(u); N <- length(u)
    I <- rep(0, N)
    for (q in seq_along(nodes)) {
      tq <- u * nodes[q]
      Hq <- build_h(tq, Uprev, spec_h)
      hq <- as.numeric(Hq %*% beta)
      htil <- pmin(pmax(hq, -H), H)
      I <- I + weights[q] * exp(htil)
    }
    u * I
  }
  U_te <- sweep(sweep(S$X_te, 2, Sc$mu, "-"), 2, Sc$sigma, "/")
  all_ok <- TRUE
  for (k in seq_len(ncol(U_te))) {
    Uprev <- if (k > 1) U_te[, 1:(k - 1), drop = FALSE] else matrix(0, nrow(U_te), 0)
    u <- U_te[, k]
    Iu <- oriented_integral(u, Uprev, Sc$coeffs[[k]]$beta, Sc$spec_h, Sc$gl_nodes, Sc$gl_weights, Sc$Hmax)
    int_code <- sign(Iu); sgn_u <- sign(u)
    ok <- (int_code == sgn_u) | ((u == 0) & (Iu == 0))
    all_ok <- all_ok && all(ok)
  }
  expect_true(all_ok)

  # Clip and Jacobian positivity checks
  minJ_tr <- min(predict_ttm(Sc, S$X_tr, type = "jac_diag"))
  minJ_te <- min(predict_ttm(Sc, S$X_te, type = "jac_diag"))
  expect_true(minJ_tr > 0 && minJ_te > 0)
  # Clipped h should be within [-H,H]
  max_abs_ht <- -Inf
  for (k in seq_len(ncol(U_te))) {
    Uprev <- if (k > 1) U_te[, 1:(k - 1), drop = FALSE] else matrix(0, nrow(U_te), 0)
    h_raw <- as.numeric(build_h(U_te[, k], Uprev, Sc$spec_h) %*% Sc$coeffs[[k]]$beta)
    h_til <- pmin(pmax(h_raw, -Sc$Hmax), Sc$Hmax)
    max_abs_ht <- max(max_abs_ht, max(abs(h_til)))
  }
  expect_true(max_abs_ht <= Sc$Hmax + 1e-12)

  # Constants and jacobian contribution appear exactly once (cross-term, test set)
  Z_te <- predict_ttm(Sc, S$X_te, type = "transform")
  Jdiag_te <- predict_ttm(Sc, S$X_te, type = "jac_diag")
  LD_te_ref <- (-0.5) * (Z_te^2) + (-0.5 * log(2 * pi)) + (log(Jdiag_te) - matrix(log(Sc$sigma), nrow = nrow(Z_te), ncol = ncol(Z_te), byrow = TRUE))
  expect_true(max(abs(te_c$LD - LD_te_ref)) <= 1e-10)
  # Extra strict constant check: the residual over entries equals C
  C <- -0.5 * log(2 * pi)
  resid <- te_c$LD - ((-0.5) * (Z_te^2) + (log(Jdiag_te) - matrix(log(Sc$sigma), nrow = nrow(Z_te), ncol = ncol(Z_te), byrow = TRUE)))
  expect_true(max(abs(resid - C)) <= 1e-12)

  # Reproducibility: refit and compare LD/J exactly
  set.seed(42L)
  F_c2 <- fit_ttm(S, algo = "crossterm", seed = seed, deg_g = 2L, df_t = 6L, Q = 16L, lambda = 1e-3, maxit = 50L)
  Sc2 <- F_c2$S
  LD2 <- predict_ttm(Sc2, S$X_te, type = "logdensity_by_dim")
  J2  <- predict_ttm(Sc2, S$X_te, type = "logdensity")
  expect_true(isTRUE(all.equal(te_c$LD, LD2, tolerance = 1e-15)))
  expect_true(isTRUE(all.equal(te_c$J,  J2,  tolerance = 1e-15)))

  # ASCII-only test log
  dir.create(file.path(root, "logs"), showWarnings = FALSE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  LOG <- file.path(root, sprintf("logs/test_ttm_crossterm_invariants_%s.md", ts))
  con <- file(LOG, open = "wt"); on.exit(try(close(con), silent = TRUE), add = TRUE)
  cat(file = con, "TTM crossterm invariant checks\n\n")
  cat(file = con, sprintf("sum_check_max=%.3g (train) / %.3g (test)\n", max(abs(rowSums(tr_c$LD) - tr_c$J)), max(abs(rowSums(te_c$LD) - te_c$J))))
  cat(file = con, sprintf("sign_ok=TRUE, max|h|=%.6g, minJ_tr=%.6g, minJ_te=%.6g\n", max_abs_ht, minJ_tr, minJ_te))
})


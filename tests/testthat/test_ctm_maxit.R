context("Cross-term maxit control and lean logging")

test_that("Higher cross.maxit does not worsen test NLL (±1e-3)", {
  root_path <- testthat::test_path("..")
  source(file.path(root_path, "02_split.R"))
  source(file.path(root_path, "R/ttm_bases.R"))
  source(file.path(root_path, "R/ttm_core.R"))
  source(file.path(root_path, "R/ttm_marginal.R"))
  source(file.path(root_path, "R/ttm_crossterm.R"))

  set.seed(777)
  N <- 60L; K <- 2L
  X <- matrix(rnorm(N * K), ncol = K)
  S <- split_data(X, seed = 777)

  # Light basis/quad for runtime
  deg_g <- 1L; df_t <- 4L; Q <- 8L; lambda <- 1e-3

  old_opt <- getOption("cross.maxit", NULL)
  on.exit({ if (is.null(old_opt)) options(cross.maxit = NULL) else options(cross.maxit = old_opt) }, add = TRUE)
  options(cross.maxit = 50L)
  fit50 <- fit_ttm(S, algo = "crossterm", deg_g = deg_g, df_t = df_t, Q = Q, lambda = lambda, seed = 777)
  options(cross.maxit = 200L)
  fit200 <- fit_ttm(S, algo = "crossterm", deg_g = deg_g, df_t = df_t, Q = Q, lambda = lambda, seed = 777)

  nll50 <- fit50$NLL_test
  nll200 <- fit200$NLL_test
  expect_lte(nll200, nll50 + 1e-3)
})

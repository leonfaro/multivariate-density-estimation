context("True half-moon log-density numerical sanity")

test_that("True half-moon density returns finite values and correct shapes", {
  source(file.path(root_path, "scripts/halfmoon_data.R"))
  source(file.path(root_path, "scripts/true_halfmoon_density.R"))
  S <- make_halfmoon_splits(n_train = 100L, n_test = 60L, noise = 0.15, seed = 5, val_frac = 0.2)
  G <- rbind(S$X_tr[1:10, , drop = FALSE], S$X_te[1:10, , drop = FALSE])
  res <- true_logdensity(G, S, Q = 16L)
  expect_true(is.numeric(res$joint))
  expect_equal(length(res$joint), nrow(G))
  expect_true(all(is.finite(res$joint)))
  expect_true(is.matrix(res$by_dim))
  expect_equal(dim(res$by_dim), c(nrow(G), 2L))
  expect_true(all(is.finite(res$by_dim)))
  expect_lte(max(abs(rowSums(res$by_dim) - res$joint)), 1e-10)
})

test_that("Conditional true density p(x|y) has sane invariants", {
  source(file.path(root_path, "scripts/halfmoon_data.R"))
  source(file.path(root_path, "scripts/true_halfmoon_density.R"))
  S <- make_halfmoon_splits(n_train = 80L, n_test = 40L, noise = 0.2, seed = 11, val_frac = 0.2)
  G <- S$X_te
  yc <- S$y_te
  res <- true_logdensity_conditional(G, yc, S, Q = 16L)
  expect_true(is.numeric(res$joint))
  expect_equal(length(res$joint), nrow(G))
  expect_true(all(is.finite(res$joint)))
  expect_true(is.matrix(res$by_dim))
  expect_equal(dim(res$by_dim), c(nrow(G), 2L))
  expect_true(all(is.finite(res$by_dim)))
  expect_lte(max(abs(rowSums(res$by_dim) - res$joint)), 1e-10)
})

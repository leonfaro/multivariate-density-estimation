library(testthat)
source('../../scripts/halfmoon_data.R')
source('../../models/ttm_cross_term.R')

set.seed(7)
S <- make_halfmoon_splits(100, 100, 0.15, 7, 0.2)
fit <- trainCrossTermMap(S, seed = 7)
M <- fit$S

test_that('Predict API shapes and sums', {
  X <- S$X_te
  LD <- predict(M, X, 'logdensity_by_dim')
  expect_equal(dim(LD), c(nrow(X), ncol(X)))
  LJ <- as.numeric(predict(M, X, 'logdensity'))
  expect_true(max(abs(rowSums(LD) - LJ)) <= 1e-10)
  expect_true(all(is.finite(LD)))
})

test_that('Eq.(22) numeric derivative matches exp(h)', {
  old <- getOption('mde.check_eq22')
  options(mde.check_eq22 = TRUE)
  on.exit(options(mde.check_eq22 = old), add = TRUE)
  # will error if violation; otherwise returns LD
  LD <- predict(M, S$X_te[1:10, ], 'logdensity_by_dim')
  expect_true(is.matrix(LD))
})


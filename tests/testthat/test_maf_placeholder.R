test_that("MAF placeholder throws", {
  S <- list(X_tr = matrix(rnorm(100*4), 100, 4))
  config <- list()
  expect_error(fit_MAF(S, config))
})

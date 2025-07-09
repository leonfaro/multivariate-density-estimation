test_that("MAF fits and predicts", {
  set.seed(1)
  S <- list(X_tr = matrix(rnorm(100*4), 100, 4),
            X_te = matrix(rnorm(100*4), 100, 4))
  config <- list()
  model <- fit_MAF(S, config, n_epoch = 1)
  lp    <- predict(model, S$X_te, "logdensity")
  expect_equal(length(lp), 100)
  expect_true(all(is.finite(lp)))
})

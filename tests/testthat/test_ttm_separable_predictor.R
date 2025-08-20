source("helper_config.R")
source("../../04_evaluation.R")

set.seed(13)
n <- 50
prep <- prepare_data(n, config)
fit <- trainSeparableMap(prep$S)

X_te <- prep$S$X_te
ld_by_dim <- predict(fit$S, X_te, "logdensity_by_dim")
ld_total <- predict(fit$S, X_te, "logdensity")

N <- nrow(X_te)
K <- ncol(X_te)

test_that("logdensity_by_dim liefert korrekte Form", {
  expect_equal(dim(ld_by_dim), c(N, K))
})

test_that("Keine NA oder Inf in den Logdichten", {
  expect_true(all(is.finite(ld_by_dim)))
  expect_true(all(is.finite(ld_total)))
})

test_that("Summen stimmen mit Gesamtlogdichte \u00fcberein", {
  expect_equal(rowSums(ld_by_dim), ld_total, tolerance = 1e-12)
})

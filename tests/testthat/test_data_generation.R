source("../../01_data_generation.R")

# Basic functionality of gen_samples

test_that("gen_samples draws matrix of correct dimension", {
  G <- setup_global()
  X <- gen_samples(G)
  expect_true(is.matrix(X))
  expect_equal(dim(X), c(G$N, length(G$config)))
  expect_false(anyNA(X))
})


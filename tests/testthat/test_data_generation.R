source("../../00_globals.R")
source("../../01_data_generation.R")

# Basic functionality of gen_samples

test_that("gen_samples draws matrix of correct dimension", {
  G <- setup_global()
  X <- gen_samples(G)
  expect_true(is.matrix(X))
  expect_equal(dim(X), c(G$n, length(G$config)))
  expect_false(anyNA(X))
})

test_that("Generate_iid_from_config returns parameter info", {
  set.seed(1)
  res <- Generate_iid_from_config(5, config, return_params = TRUE)
  expect_true(is.matrix(res$X))
  expect_type(res$params, "list")
  expect_equal(names(res$params[[2]]), "rate")
  expect_equal(names(res$params[[3]]), c("shape1", "shape2"))
  expect_equal(names(res$params[[4]]), c("shape", "scale"))
  expect_false(anyNA(unlist(res$params[[2]])))
})


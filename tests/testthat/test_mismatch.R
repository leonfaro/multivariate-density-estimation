library(testthat)

## run tests from repository root
if (file.exists('../../run5.R')) {
  capture.output(source('../../run5.R', chdir = TRUE))
} else {
  skip('run5.R not available')
}

test_that('joint mismatches are finite', {
  expect_true(is.finite(forest_mismatch))
  expect_true(is.finite(kernel_mismatch))
  expect_true(is.finite(copula_mismatch))
  expect_true(all(is.finite(eval_tab$delta_param)))
  expect_true(all(is.finite(eval_tab$delta_forest)))
  expect_true(all(is.finite(eval_tab$delta_kernel)))
  expect_true(all(is.finite(eval_tab$delta_dvine)))
})

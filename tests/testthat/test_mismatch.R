library(testthat)

## run tests from repository root
capture.output(source('../../run5.R', chdir = TRUE))

test_that('joint mismatches are finite', {
  expect_true(is.finite(forest_mismatch))
  expect_true(is.finite(kernel_mismatch))
  expect_true(is.finite(copula_mismatch))
  expect_true(all(is.finite(eval_tab$delta_param)))
  expect_true(all(is.finite(eval_tab$delta_forest)))
  expect_true(all(is.finite(eval_tab$delta_kernel)))
  expect_true(all(is.finite(eval_tab$delta_dvine)))
})

test_that('determinant of the Jacobian is positive', {
  jac_train <- exp(det_J(logd_train))
  jac_test  <- exp(det_J(logd_test))
  expect_true(all(jac_train > 0))
  expect_true(all(jac_test  > 0))
})

test_that('diagonal Jacobian entries are positive', {
  diag_train <- exp(logd_train)
  diag_test  <- exp(logd_test)
  expect_true(all(diag_train > 0))
  expect_true(all(diag_test  > 0))
})

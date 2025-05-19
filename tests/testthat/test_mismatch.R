library(testthat)

capture.output(source('run5.R'))

test_that('joint mismatches vanish', {
  expect_lt(abs(forest_mismatch), 1e-6)
  expect_lt(abs(kernel_mismatch), 1e-6)
  expect_lt(abs(copula_mismatch), 1e-6)
  expect_true(all(abs(eval_tab$delta_param) < 1e-6))
  expect_true(all(abs(eval_tab$delta_forest) < 1e-6))
  expect_true(all(abs(eval_tab$delta_kernel) < 1e-6))
  expect_true(all(abs(eval_tab$delta_dvine) < 1e-6))
})

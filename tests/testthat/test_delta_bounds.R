library(testthat)

if (file.exists('../../run5.R')) {
  capture.output(source('../../run5.R', chdir = TRUE))
} else {
  skip('run5.R not available')
}

test_that('delta columns bounded by 1', {
  expect_true(all(is.finite(eval_tab$delta_ll_joint)))
  expect_true(all(is.finite(eval_tab$delta_ll_trtf)))
  expect_true(all(is.finite(eval_tab$delta_ll_kernel)))
  expect_true(all(abs(eval_tab$delta_ll_joint) <= 1))
  expect_true(all(abs(eval_tab$delta_ll_trtf) <= 1))
  expect_true(all(abs(eval_tab$delta_ll_kernel) <= 2))
})

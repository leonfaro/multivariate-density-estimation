library(testthat)

if (file.exists('../../run5.R')) {
  capture.output(source('../../run5.R', chdir = TRUE))
} else {
  skip('run5.R not available')
}

test_that('delta columns finite', {
  expect_true(all(is.finite(eval_tab$delta_ll_param)))
  expect_true(all(is.finite(eval_tab$delta_ll_trtf)))
  expect_true(all(is.finite(eval_tab$delta_ll_kernel)))
})

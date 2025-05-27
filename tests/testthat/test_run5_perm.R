library(testthat)

if (file.exists('../../run5.R')) {
  capture.output(source('../../run5.R', chdir = TRUE))
} else {
  skip('run5.R not available')
}

test_that('eval_tab_nat and eval_tab_perm exist and are finite', {
  expect_true(exists('eval_tab_nat'))
  expect_true(exists('eval_tab_perm'))
  for (tbl in list(eval_tab_nat, eval_tab_perm)) {
    num_cols <- sapply(tbl, is.numeric)
    expect_true(all(is.finite(as.matrix(tbl[, num_cols]))))
  }
})

test_that("setup_global returns expected list", {
  out <- setup_global()
  expect_type(out, "list")
  expect_equal(out$n, 500)
  expect_equal(length(out$config), length(config))
  expect_equal(vapply(out$config, `[[`, "distr", FUN.VALUE = ""),
               vapply(config, `[[`, "distr", FUN.VALUE = ""))
  expect_equal(out$seed, 42)
  expect_equal(out$split_ratio, 0.5)
  expect_equal(out$p_max, 6)
  expect_equal(out$h_grid, seq_len(out$p_max))
  expect_equal(out$model_ids, "TRUE")
})

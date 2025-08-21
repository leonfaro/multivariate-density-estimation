source(file.path(testthat::test_path("..", ".."), "main.R"), chdir = TRUE)

test_that("main handles halfmoon2d dataset", {
  withr::local_envvar(c(DATASET="halfmoon2d", N_TRAIN="20", N_TEST="20", NOISE="0.1", SEED="7"))
  f <- "results/splits_halfmoon2d_seed007.rds"
  unlink(f)
  on.exit(unlink(f), add = TRUE)
  expect_output(main(), "\\[DATASET halfmoon2d\\]")
  expect_true(file.exists(f))
})

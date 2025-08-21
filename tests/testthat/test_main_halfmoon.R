source(file.path(testthat::test_path("..", ".."), "main.R"), chdir = TRUE)

test_that("main handles halfmoon2d dataset", {
  withr::local_envvar(c(DATASET="halfmoon2d", N_TRAIN="20", N_TEST="20", NOISE="0.1", SEED="7"))
  f_split <- "results/splits_halfmoon2d_seed007.rds"
  f_csv <- "results/nll_halfmoon_seed007.csv"
  unlink(c(f_split, f_csv))
  on.exit(unlink(c(f_split, f_csv)), add = TRUE)
  expect_output(main(), "\\[DATASET halfmoon2d\\]")
  expect_true(file.exists(f_split))
  expect_true(file.exists(f_csv))
  df <- read.csv(f_csv, stringsAsFactors = FALSE)
  expect_equal(df$model, c("true", "trtf", "ttm", "ttm_sep", "ttm_cross"))
  expect_true(all(is.finite(as.matrix(df[,-1]))))
})

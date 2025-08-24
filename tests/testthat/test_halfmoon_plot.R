source("../../00_globals.R")
source("../../models/true_model.R")
source("../../scripts/halfmoon_data.R")
source("../../scripts/halfmoon_plot.R")

set.seed(1)
S <- make_halfmoon_splits(20, 10, 0.1, seed = 1)
mods <- list(
  true = list(theta = list(c(0, 1), c(0, 1))),
  trtf = structure(list(), class = "dummy"),
  ttm = structure(list(), class = "dummy"),
  ttm_sep = structure(list(), class = "dummy"),
  ttm_cross = structure(list(), class = "dummy")
)
predict.dummy <- function(object, x, ...) -rowSums(x^2)
assign("predict.dummy", predict.dummy, envir = .GlobalEnv)

test_that("draw_points returns the total number of points", {
  tmp <- tempfile(fileext = ".png")
  png(tmp)
  plot(NA, xlim = c(-3, 3), ylim = c(-3, 3))
  n <- draw_points(S)
  dev.off()
  expect_equal(n, nrow(S$X_tr) + nrow(S$X_val) + nrow(S$X_te))
})

test_that("global_quantiles policy yields common levels", {
  skip_on_cran()
  tmp <- tempfile(fileext = ".png")
  png(tmp)
  res <- .draw_panels(mods, S, grid_n = 160, levels_policy = "global_quantiles")
  dev.off()
  expect_equal(res$grid_points, 160^2)
  expect_true(res$all_finite)
  expect_length(res$levels, 3)
})

test_that("per_model policy returns separate levels", {
  skip_on_cran()
  tmp <- tempfile(fileext = ".png")
  png(tmp)
  res <- .draw_panels(mods, S, grid_n = 160, levels_policy = "per_model")
  dev.off()
  expect_equal(res$grid_points, 160^2)
  expect_true(res$all_finite)
  expect_equal(length(res$levels), 5)
  expect_true(all(vapply(res$levels, length, 1L) == 3))
})

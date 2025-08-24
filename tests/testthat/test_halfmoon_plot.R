source("../../scripts/halfmoon_data.R")
source("../../scripts/halfmoon_plot.R")

set.seed(1)
S <- make_halfmoon_splits(20, 10, 0.1, seed = 1)

test_that("draw_points returns the total number of points", {
  tmp <- tempfile(fileext = ".png")
  png(tmp)
  plot(NA, xlim = c(-3, 3), ylim = c(-3, 3))
  n <- draw_points(S)
  dev.off()
  expect_equal(n, nrow(S$X_tr) + nrow(S$X_val) + nrow(S$X_te))
})

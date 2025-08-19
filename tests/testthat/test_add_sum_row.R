source("../../04_evaluation.R")

# Test that add_sum_row ignores NA values in numeric columns

set.seed(1)

test_that("add_sum_row summiert numerische Spalten trotz NA", {
  tab <- data.frame(dim = c("1", "2"), a = c(1, NA), b = c(2, 3))
  res <- add_sum_row(tab, label = "tot")
  expect_equal(res[3, "a"], 1)
  expect_equal(res[3, "b"], 5)
  expect_equal(res[3, "dim"], "tot")
})

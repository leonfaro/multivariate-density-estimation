source("../../replicate_code.R")
results_table <- matrix(1:4, nrow = 2)
replicate_code_scripts("../../main.R", "tmp_repl.txt", env = environment())
out <- readLines("tmp_repl.txt")

test_that("replicate_code captures main and table", {
expect_true(any(grepl("Begin .*main.R", out)))
  expect_true(any(grepl("Final results table", out)))
  tbl <- tail(out, n = 2)
  nums <- as.numeric(unlist(strsplit(gsub("[^0-9]", " ", paste(tbl, collapse=" ")), " ")))
  nums <- nums[!is.na(nums)]
  expect_true(all(is.finite(nums)))
})

unlink("tmp_repl.txt")


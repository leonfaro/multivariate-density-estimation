source("helper_config.R")
source("../../EDA.R")
source("../../models/trtf_model.R")
source("../../models/ks_model.R")
source("../../models/true_model.R")
source("../../02_split.R")

set.seed(6)
res <- run_pipeline(20, config, c(3,2,1,4))

test_that("run_pipeline returns comprehensive list", {
  expect_true(all(c("data","models","tables","scatter_data","plots","times","runtime") %in% names(res)))
  expect_length(res$plots$scatter, 4)
  expect_true(all(vapply(res$plots$scatter, inherits, logical(1), "recordedplot")))
  expect_s3_class(res$tables$combined, "knitr_kable")
})

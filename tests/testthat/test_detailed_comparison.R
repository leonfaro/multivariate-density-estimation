library(testthat)
source("../../00_setup.R", chdir = TRUE)

test_that("save_detailed_comparison_data outputs CSV", {
  set.seed(42)
  data <- generate_data(N_train = 10, N_test = 10)
  X <- data$train$sample$X_pi
  param_ests <- lapply(seq_len(K), function(k) {
    family_spec <- dist_registry[[config[[k]]$distr]]
    n_par <- length(family_spec$param_names)
    rep(0, n_par * if (k == 1) 1 else (k - 1) + 1)
  })
  tmp <- tempfile()
  dir.create(tmp)
  save_detailed_comparison_data(
    data_matrix = X,
    param_ests_list = param_ests,
    config_list = config,
    dist_registry_obj = dist_registry,
    output_dir_path = tmp,
    data_label_str = "train"
  )
  f1 <- file.path(tmp, paste0("dim1_", config[[1]]$distr, "_params_logpdf_train.csv"))
  expect_true(file.exists(f1))
  df <- read.csv(f1)
  expect_equal(nrow(df), nrow(X))
  expect_false(any(is.na(df)))
})

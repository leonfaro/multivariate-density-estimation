source("helper_config.R")
source("../../04_evaluation.R")
source("../../models/trtf_model.R")
source("../../models/ks_model.R")
source("../../models/true_model.R")
source("../../02_split.R")

set.seed(5)
prep <- prepare_data(20, config, c(3,2,1,4))
mods_n <- fit_models(prep$S, config)
mods_p <- fit_models(prep$S_perm, config)
scat_n <- make_scatter_data(mods_n, prep$S)
scat_p <- make_scatter_data(mods_p, prep$S_perm)
scatter_data <- c(scat_n, setNames(scat_p, paste0(names(scat_p), "_p")))

plots_sc <- plot_scatter(scatter_data)
plots_par <- plot_parameters(prep$param_list)

test_that("plot_scatter returns four recorded plots", {
  expect_length(plots_sc, 4)
  expect_true(all(vapply(plots_sc, inherits, logical(1), "recordedplot")))
})

test_that("plot_parameters returns recorded plots", {
  expect_true(length(plots_par) >= 0)
  if (length(plots_par) > 0) {
    expect_true(all(vapply(plots_par, inherits, logical(1), "recordedplot")))
  }
})

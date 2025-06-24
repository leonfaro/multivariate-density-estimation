source("helper_config.R")
source("../../04_evaluation.R")
source("../../models/trtf_model.R")
source("../../models/ks_model.R")
source("../../models/true_model.R")
source("../../02_split.R")
sys.source("../../main.R", envir = environment(), chdir = TRUE)

set.seed(7)
prep <- prepare_data(20, config, c(3,2,1,4))
mods_n <- fit_models(prep$S, config)
mods_p <- fit_models(prep$S_perm, config)
scat_n <- make_scatter_data(mods_n, prep$S)
scat_p <- make_scatter_data(mods_p, prep$S_perm)
scatter_data <- c(scat_n, setNames(scat_p, paste0(names(scat_p), "_p")))

plt <- plot_scatter_matrix(scatter_data)

test_that("plot_scatter_matrix returns recorded plot", {
  expect_s3_class(plt, "recordedplot")
})

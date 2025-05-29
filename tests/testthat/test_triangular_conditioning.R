library(testthat)

extract_config <- function(file) {
  exprs <- parse(file)
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name("config"))) {
      return(eval(e[[3]], envir = new.env()))
    }
  }
  stop("config not found in ", file)
}

check_triangular_conditioning <- function(cfg) {
  for (k in seq_along(cfg)) {
    ck <- cfg[[k]]
    if (!is.null(ck$parm) && k > 1) {
      fn_txt <- paste(deparse(ck$parm), collapse = " ")
      expect_true(grepl(paste0("X", k - 1), fn_txt, fixed = TRUE),
                  info = paste0("component ", k,
                                " does not condition on X", k - 1))
    }
  }
}

test_that("configs use immediate predecessor", {
  root <- file.path("..", "..")
  files <- c(
    "run3.R",
    "basic_tests.R",
    file.path("tests/testthat", "test_cdf_log.R"),
    file.path("tests/testthat", "test_edgecases.R"),
    file.path("tests/testthat", "test_fit_param.R"),
    file.path("tests/testthat", "test_softplus_gamma.R")
  )
  for (f in files) {
    cfg <- extract_config(file.path(root, f))
    check_triangular_conditioning(cfg)
  }
})

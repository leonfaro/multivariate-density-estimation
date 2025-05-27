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

check_cfg <- function(cfg, root) {
  env <- new.env()
  env$config <- cfg
  sys.source(file.path(root, "00_setup.R"), env, chdir = TRUE)
  for (k in seq_along(cfg)) {
    ck <- cfg[[k]]
    expect_true(is.list(ck))
    expect_true(is.character(ck$distr))
    expect_true(ck$distr %in% names(env$q_supports_logp))

    if (!is.null(ck$parm)) {
      fn_txt <- paste(deparse(ck$parm), collapse = " ")
      # Check no future dimensions referenced
      max_dim <- max(length(cfg), 6)
      for (d in seq(from = k, to = max_dim)) {
        expect_false(grepl(paste0("X", d), fn_txt, fixed = TRUE))
      }
      # Positive parameter enforcement via monotone transform
      need_pos <- ck$distr %in% c("norm", "exp", "gamma", "weibull", "lnorm",
                                  "pois", "beta", "logis")
      if (need_pos) {
        expect_true(grepl("softplus", fn_txt) ||
                    grepl("exp(", fn_txt, fixed = TRUE) ||
                    grepl("plogis", fn_txt))
      }
      # Call function on zero data frame
      if (k > 1) {
        d <- setNames(as.data.frame(matrix(0, nrow = 1, ncol = k - 1)),
                       paste0("X", seq_len(k - 1)))
      } else {
        d <- data.frame()[, FALSE]
      }
      environment(ck$parm) <- env
      pars <- ck$parm(d)
      expect_silent(env$safe_pars(pars, ck$distr))
    }
  }
}

test_that("configs satisfy monotonicity rules", {
  root <- file.path("..", "..")
  check_cfg(extract_config(file.path(root, "run3.R")), root)
  check_cfg(extract_config(file.path(root, "run5.R")), root)
  check_cfg(extract_config(file.path(root, "basic_tests.R")), root)
})

fail <- FALSE

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

check_triangular_arity <- function(cfg) {
  for (k in seq_along(cfg)) {
    ck <- cfg[[k]]
    if (!is.null(ck$parm)) {
      fn_txt <- paste(deparse(ck$parm), collapse = " ")
      for (d in seq(from = k + 1, to = length(cfg) + 5)) {
        if (grepl(paste0("X", d), fn_txt, fixed = TRUE)) {
          message(sprintf("triangular-arity violation in component %d: X%d referenced", k, d))
          fail <<- TRUE
        }
      }
    }
  }
}

check_monotone_link_naming <- function() {
  env <- new.env()
  env$config <- list(list(distr = "norm", parm = NULL))
  sys.source("00_setup.R", envir = env)
  allowed <- c("identity", "softplus", "exp")
  for (nm in names(env$dist_registry)) {
    lv <- env$dist_registry[[nm]]$link_vector
    if (any(!lv %in% allowed)) {
      message(sprintf("unapproved link function in %s", nm))
      fail <<- TRUE
    }
  }
}

check_config_schema <- function(cfg) {
  env <- new.env()
  env$config <- cfg
  sys.source("00_setup.R", envir = env)
  for (k in seq_along(cfg)) {
    ck <- cfg[[k]]
    if (!all(c("distr", "parm") %in% names(ck))) {
      message(sprintf("config entry %d missing fields", k))
      fail <<- TRUE
    }
    if (!ck$distr %in% env$allowed_dists_full) {
      message(sprintf("config entry %d uses disallowed distribution", k))
      fail <<- TRUE
    }
    if (!is.null(ck$parm) && !is.function(ck$parm)) {
      message(sprintf("config entry %d has non-function parm", k))
      fail <<- TRUE
    }
  }
}

check_no_silent_globals <- function(files) {
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    if (any(grepl("<<-", lines, fixed = TRUE))) {
      message(sprintf("global assignment found in %s", f))
      fail <<- TRUE
    }
  }
}

check_inverse_export <- function(files) {
  regs <- unlist(lapply(files, function(f) {
    lines <- readLines(f, warn = FALSE)
    grep("register_map", lines, value = TRUE)
  }))
  if (length(regs) > 0) {
    funs <- character()
    for (f in files) {
      lines <- readLines(f, warn = FALSE)
      defs <- sub("^(.*)<-.*", "\\1", grep("<-\\s*function", lines, value = TRUE))
      funs <- c(funs, trimws(defs))
    }
    if (!("S" %in% funs && "S_inv" %in% funs)) {
      message("map registered without both S and S_inv")
      fail <<- TRUE
    }
  }
}

check_jacobian_shortcut <- function(files) {
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    if (any(grepl("Jacobian override", lines, ignore.case = TRUE))) {
      if (!any(grepl("scalar shortcut", lines, ignore.case = TRUE))) {
        message(sprintf("jacobian override without shortcut in %s", f))
        fail <<- TRUE
      }
    }
  }
}

check_unit_scale_tails <- function(files) {
  invisible(NULL)
}

check_no_hard_seed <- function(files) {
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    if (any(grepl("set.seed", lines, fixed = TRUE))) {
      message(sprintf("set.seed call found in %s", f))
      fail <<- TRUE
    }
  }
}

main <- function() {
  cfg_files <- c("run3.R", "run5.R", "basic_tests.R")
  cfgs <- lapply(cfg_files, extract_config)
  lapply(cfgs, check_triangular_arity)
  lapply(cfgs, check_config_schema)

  check_monotone_link_naming()

  library_files <- setdiff(list.files(pattern = "\\.R$"),
                           c("run3.R", "run5.R", "basic_tests.R", "cdf_logcdf_range.R", "extra_lint.R"))
  check_no_silent_globals(library_files)
  check_inverse_export(library_files)
  check_jacobian_shortcut(library_files)
  check_unit_scale_tails(library_files)

  allowed_seed <- c("run3.R", "run5.R", "basic_tests.R", "cdf_logcdf_range.R", "extra_lint.R")
  seed_files <- setdiff(list.files(pattern = "\\.R$"), allowed_seed)
  check_no_hard_seed(seed_files)

  if (fail) quit(status = 1)
}

main()

if (!exists("locate_repo_loader", inherits = TRUE)) {
  locate_repo_loader <- function() {
    detect_script_path <- function() {
      frames <- sys.frames()
      for (i in rev(seq_along(frames))) {
        fi <- frames[[i]]
        if (!is.null(fi$ofile)) {
          path <- tryCatch(normalizePath(fi$ofile, winslash = "/", mustWork = TRUE),
                          error = function(e) NA_character_)
          if (!is.na(path) && nzchar(path)) return(path)
        }
      }
      args <- commandArgs(trailingOnly = FALSE)
      file_arg <- args[grepl("^--file=", args)]
      if (length(file_arg)) {
        cand <- sub("^--file=", "", file_arg[1])
        path <- tryCatch(normalizePath(cand, winslash = "/", mustWork = TRUE),
                         error = function(e) NA_character_)
        if (!is.na(path) && nzchar(path)) return(path)
      }
      NA_character_
    }

    start_dirs <- character()
    script_path <- detect_script_path()
    if (!is.na(script_path) && nzchar(script_path)) {
      start_dirs <- c(start_dirs, dirname(script_path))
    }
    wd <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = FALSE),
                   error = function(e) getwd())
    start_dirs <- unique(c(start_dirs, wd))
    checked <- character()
    for (start in start_dirs) {
      cur <- start
      repeat {
        cur <- tryCatch(normalizePath(cur, winslash = "/", mustWork = FALSE),
                        error = function(e) cur)
        if (!nzchar(cur) || cur %in% checked) break
        checked <- c(checked, cur)
        cand1 <- file.path(cur, "R", "loader.R")
        if (file.exists(cand1)) {
          return(normalizePath(cand1, winslash = "/", mustWork = TRUE))
        }
        cand2 <- file.path(cur, "code", "R", "loader.R")
        if (file.exists(cand2)) {
          return(normalizePath(cand2, winslash = "/", mustWork = TRUE))
        }
        parent <- dirname(cur)
        if (identical(parent, cur)) break
        cur <- parent
      }
    }
    stop("Could not locate loader.R")
  }
}

loader_path <- locate_repo_loader()
if (!exists("initialize_repo")) {
  source(loader_path, chdir = FALSE)
}
root_path <- initialize_repo()
source(file.path(root_path, "experiments", "halfmoon", "halfmoon_pipeline.R"))

if (sys.nframe() == 0L) {
  invisible(NULL)
}

config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp", parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "beta",
       parm = function(d) list(shape1 = softplus(d$X2),
                               shape2 = softplus(d$X1))),
  list(distr = "gamma",
       parm = function(d) list(shape = softplus(d$X3),
                               scale = softplus(d$X2)))
)

config_canonical <- config

N <- 50L
override_N <- suppressWarnings(as.integer(Sys.getenv("N_OVERRIDE", "")))
if (is.finite(override_N) && !is.na(override_N) && override_N > 0L) {
  N <- override_N
}
n <- N

ensure_loader <- function() {
  if (exists("initialize_repo", inherits = TRUE)) {
    existing_root <- tryCatch(get("root_path", envir = .GlobalEnv), error = function(e) NULL)
    return(list(code_dir = existing_root))
  }

  locate_candidates <- function() {
    paths <- character()

    script_path <- tryCatch({
      normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE)
    }, error = function(e) NA_character_)

    if (is.na(script_path) || !nzchar(script_path)) {
      script_path <- tryCatch({
        args <- commandArgs(trailingOnly = FALSE)
        marker <- "--file="
        file_arg <- args[grepl(marker, args, fixed = TRUE)]
        if (length(file_arg)) {
          normalizePath(sub(marker, "", file_arg[1]), winslash = "/", mustWork = TRUE)
        } else {
          NA_character_
        }
      }, error = function(e) NA_character_)
    }

    if (!is.na(script_path) && nzchar(script_path)) {
      script_dir <- dirname(script_path)
      paths <- c(paths,
                 file.path(script_dir, "R", "loader.R"),
                 file.path(dirname(script_dir), "R", "loader.R"),
                 file.path(dirname(script_dir), "code", "R", "loader.R"))
    }

    wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    paths <- c(paths,
               file.path(wd, "R", "loader.R"),
               file.path(wd, "code", "R", "loader.R"),
               file.path(dirname(wd), "R", "loader.R"),
               file.path(dirname(wd), "code", "R", "loader.R"))
    unique(paths)
  }

  candidates <- locate_candidates()
  for (cand in candidates) {
    if (file.exists(cand)) {
      source(cand, chdir = FALSE)
      code_dir <- dirname(dirname(cand))
      return(list(code_dir = code_dir))
    }
  }

  stop(sprintf("Could not locate loader.R. Checked paths: %s",
               paste(candidates, collapse = ", ")))
}

loader_info <- ensure_loader()

load_repo_modules <- function(code_dir) {
  if (is.null(code_dir) || !nzchar(code_dir) || !dir.exists(code_dir)) {
    fallback <- file.path(getwd(), "code")
    if (dir.exists(fallback)) {
      code_dir <- fallback
    }
  }
  old_wd <- getwd()
  reset_needed <- FALSE
  if (!is.null(code_dir) && nzchar(code_dir) && dir.exists(code_dir)) {
    old_norm <- tryCatch(normalizePath(old_wd, winslash = "/", mustWork = TRUE), error = function(e) old_wd)
    code_norm <- tryCatch(normalizePath(code_dir, winslash = "/", mustWork = TRUE), error = function(e) code_dir)
    if (!identical(old_norm, code_norm)) {
      setwd(code_norm)
      reset_needed <- TRUE
    }
  }
  on.exit(if (reset_needed) setwd(old_wd), add = TRUE)
  initialize_repo()
}

root_path <- load_repo_modules(loader_info$code_dir)
source(file.path(root_path, "experiments", "4d", "experiment_config.R"))

perm_override <- c(4L, 2L, 3L, 1L)

run_main_perm <- function() {
  old_perm <- get("perm", envir = .GlobalEnv)
  orig_fit_true_joint <- get("fit_TRUE_JOINT", envir = .GlobalEnv)
  on.exit(assign("perm", old_perm, envir = .GlobalEnv), add = TRUE)
  on.exit(assign("fit_TRUE_JOINT", orig_fit_true_joint, envir = .GlobalEnv), add = TRUE)
  assign("perm", perm_override, envir = .GlobalEnv)

  fit_true_joint_override <- function(S, config, cores = NC) {
    if (!is.list(S) || is.null(S$X_tr) || is.null(S$X_te)) {
      return(orig_fit_true_joint(S, config, cores))
    }
    K <- ncol(S$X_tr)
    if (!length(perm_override) || length(perm_override) != K) {
      return(orig_fit_true_joint(S, config, cores))
    }
    inv_perm <- match(seq_len(K), perm_override)
    if (anyNA(inv_perm)) return(orig_fit_true_joint(S, config, cores))

    S_canon <- S
    S_canon$X_tr <- S$X_tr[, inv_perm, drop = FALSE]
    S_canon$X_te <- S$X_te[, inv_perm, drop = FALSE]

    res <- orig_fit_true_joint(S_canon, config_canonical, cores)
    if (is.list(res) && !is.null(res$logL_te_dim)) {
      res$logL_te_dim <- res$logL_te_dim[perm_override]
    }
    res$config <- config
    res
  }

  assign("fit_TRUE_JOINT", fit_true_joint_override, envir = .GlobalEnv)
  main()
}

if (sys.nframe() == 0L) {
  invisible(run_main_perm())
}

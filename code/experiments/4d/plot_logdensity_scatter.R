#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Locate this script and work from its directory so relative paths match
locate_script_dir <- function() {
  frame_files <- vapply(sys.frames(), function(fr) {
    if (!is.null(fr$ofile)) return(fr$ofile)
    NA_character_
  }, character(1), USE.NAMES = FALSE)
  cand <- frame_files[!is.na(frame_files)][1]
  if (!is.na(cand)) {
    return(dirname(normalizePath(cand, winslash = "/", mustWork = TRUE)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  marker <- "--file="
  file_arg <- args[grepl(marker, args, fixed = TRUE)]
  if (length(file_arg)) {
    return(dirname(normalizePath(sub(marker, "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- locate_script_dir()
code_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
old_wd <- getwd()
setwd(code_root)
on.exit(setwd(old_wd), add = TRUE)

Sys.setenv(N_OVERRIDE = "250")

source(file.path(script_dir, "main.R"), local = FALSE)

if (!exists("perm", envir = .GlobalEnv) || !is.numeric(get("perm", envir = .GlobalEnv))) {
  assign("perm", seq_len(length(config)), envir = .GlobalEnv)
}

seed <- as.integer(Sys.getenv("SEED", "42"))
set.seed(seed)

prep <- prepare_data(n, config, seed = seed)
S0 <- prep$S
K <- ncol(S0$X_tr)
perm_use <- get("perm", envir = .GlobalEnv)
if (!is.numeric(perm_use) || length(perm_use) != K) {
  perm_use <- seq_len(K)
  assign("perm", perm_use, envir = .GlobalEnv)
}
S <- list(
  X_tr = S0$X_tr[, perm_use, drop = FALSE],
  X_te = S0$X_te[, perm_use, drop = FALSE]
)
cfg <- config[perm_use]

fit_safe <- function(name, expr) {
  tryCatch(expr, error = function(e) {
    message(sprintf("[WARN] %s unavailable: %s", name, e$message))
    NULL
  })
}

mod_true <- fit_TRUE(S, cfg)
mod_true_joint <- fit_safe("True (Joint)", fit_TRUE_JOINT(S, cfg))
mod_trtf <- fit_safe("Random Forest", fit_TRTF(S, cfg, seed = seed))
mod_cop <- fit_safe("Copula NP", fit_copula_np(S, seed = seed))
mod_ttm <- fit_safe("Marginal Map", fit_ttm(S, algo = "marginal", seed = seed))
mod_ttm_sep <- fit_safe("Separable Map", fit_ttm(S, algo = "separable", seed = seed))

n_te <- nrow(S$X_te)
logd_true_joint <- fit_safe("true_joint_logdensity_by_dim", true_joint_logdensity_by_dim(cfg, S$X_te))
if (is.null(logd_true_joint) || !is.matrix(logd_true_joint)) {
  stop("Failed to compute true joint log-densities; cannot plot.")
}

logd_true_marg <- matrix(NA_real_, nrow = n_te, ncol = K)
for (k in seq_len(K)) {
  logd_true_marg[, k] <- .log_density_vec(S$X_te[, k], cfg[[k]]$distr, mod_true$theta[[k]])
}

predict_safe <- function(name, expr) {
  tryCatch(expr, error = function(e) {
    message(sprintf("[WARN] %s prediction failed: %s", name, e$message))
    NULL
  })
}

logd_trtf <- if (!is.null(mod_trtf)) predict_safe("Random Forest", predict(mod_trtf, S$X_te, type = "logdensity_by_dim")) else NULL

logd_ttm <- NULL
if (!is.null(mod_ttm)) {
  logd_ttm <- predict_safe("Marginal Map", predict_ttm(mod_ttm$S, S$X_te, type = "logdensity_by_dim"))
  if (is.null(logd_ttm)) {
    logd_ttm <- predict_safe("Marginal Map", predict(mod_ttm$S, S$X_te, type = "logdensity_by_dim"))
  }
}

logd_ttm_sep <- NULL
if (!is.null(mod_ttm_sep)) {
  logd_ttm_sep <- predict_safe("Separable Map", predict_ttm(mod_ttm_sep$S, S$X_te, type = "logdensity_by_dim"))
  if (is.null(logd_ttm_sep)) {
    logd_ttm_sep <- predict_safe("Separable Map", predict(mod_ttm_sep$S, S$X_te, type = "logdensity_by_dim"))
  }
}

logd_cop <- if (!is.null(mod_cop)) predict_safe("Copula NP", predict(mod_cop, S$X_te, type = "logdensity_by_dim")) else NULL

model_mats <- list(
  "True (marginal)" = logd_true_marg,
  "True (Joint)" = logd_true_joint,
  "Random Forest" = logd_trtf,
  "Marginal Map" = logd_ttm,
  "Separable Map" = logd_ttm_sep,
  "Copula NP" = logd_cop
)

model_order <- names(model_mats)
valid_models <- list()
for (nm in model_order) {
  mat <- model_mats[[nm]]
  if (is.matrix(mat) && all(dim(mat) == c(n_te, K))) {
    valid_models[[nm]] <- mat
  }
}
if (!length(valid_models)) {
  stop("No valid model predictions available for plotting.")
}

joint_true <- rowSums(logd_true_joint)
model_joint <- lapply(valid_models, rowSums)

plot_dir <- file.path(script_dir, "results")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
plot_file <- file.path(plot_dir, "logdensity_joint_N250.png")

png(plot_file, width = 8, height = 7, units = "in", res = 600)
par(mfrow = c(3, 2), mar = c(2, 2, 2, 1), oma = c(4, 4, 2, 1))

point_col <- grDevices::rgb(0, 0, 0, alpha = 0.35)
line_col <- "#C00000"
sm_col <- "#0052A5"
lims <- c(-10, 5)

for (idx in seq_along(valid_models)) {
  model_name <- names(valid_models)[idx]
  vals <- model_joint[[idx]]
  ok <- is.finite(joint_true) & is.finite(vals)
  plot(joint_true[ok], vals[ok], pch = 16, cex = 0.45,
       col = point_col, xlim = lims, ylim = lims,
       xlab = "", ylab = "")
  abline(0, 1, col = line_col, lty = 2, lwd = 1)
  if (sum(ok) > 10L) {
    sm <- stats::lowess(joint_true[ok], vals[ok], f = 0.8)
    lines(sm, col = sm_col, lwd = 1.1)
  }
}
if (length(valid_models) < 6L) {
  for (j in seq_len(6L - length(valid_models))) {
    plot.new()
  }
}

mtext("True joint log-density", side = 1, line = 2, outer = TRUE)
mtext("Model joint log-density", side = 2, line = 2.2, outer = TRUE)

dev.off()

cat(sprintf("Saved plot to %s\n", plot_file))

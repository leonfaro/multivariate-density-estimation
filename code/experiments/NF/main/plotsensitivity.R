# scripts/plot_n_sensitivity.R
# Base R script to plot N-sensitivity (NLL vs N) with error bars per dataset.
# Reads TRTF/NF-style result *.txt files from GAS/POWER/HEPMASS/MINIBOONE folders
# and extracts "sum:" as mean log-likelihood. Error bars = standard error over
# seeds for each (dataset, N). Output PNG is written to the results root.
# Usage:
#   Rscript plotsensitivity.R [RESULTS_ROOT]
# Default RESULTS_ROOT is auto-detected relative to this script.

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
default_results_root <- file.path(repo_root(dirname(dirname(loader_path))),
                                  "experiments", "NF", "main", "results")

stderr <- function(x) {
  if (!length(x)) return(NA_real_)
  if (length(x) == 1L) return(0)
  s <- stats::sd(x)
  if (is.na(s) || !is.finite(s)) return(0)
  s / sqrt(length(x))
}

extract_nll_from_file <- function(path) {
  bn <- basename(path)
  base <- sub("\\.txt$", "", bn)
  toks <- strsplit(base, "_", fixed = TRUE)[[1]]
  seed_idx <- which(grepl("^seed[0-9]+$", toks))
  n_idx <- which(grepl("^N[0-9]+$", toks))
  if (length(seed_idx)) {
    seed <- suppressWarnings(as.integer(sub("^seed", "", toks[seed_idx[1L]])))
  } else {
    seed <- NA_integer_
  }
  if (length(n_idx)) {
    N <- suppressWarnings(as.integer(sub("^N", "", toks[n_idx[1L]])))
  } else {
    N <- suppressWarnings(as.integer(sub(".*_N([0-9]+).*", "\\1", bn)))
  }
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  sum_line <- if (length(lines)) grep("^\\s*sum\\s*:", lines, value = TRUE, ignore.case = TRUE) else character(0)
  ll <- NA_real_
  if (length(sum_line)) {
    val <- sub("^.*sum\\s*:\\s*", "", sum_line[1L])
    ll <- suppressWarnings(as.numeric(val))
  }
  nll <- if (is.finite(ll)) -ll else NA_real_
  list(N = N, seed = seed, nll = nll)
}

args <- commandArgs(trailingOnly = TRUE)
results_root <- if (length(args) >= 1L && nzchar(args[1L])) args[1L] else default_results_root
results_root <- normalizePath(results_root, mustWork = FALSE)
if (!dir.exists(results_root)) {
  stop(sprintf("Results root not found: %s", results_root))
}

Ns <- c(25L, 50L, 100L, 250L, 500L, 1000L, 2500L, 5000L)
datasets <- c(gas = "Gas", power = "Power", hepmass = "Hepmass", miniboone = "Miniboone")

records <- list()
for (ds in names(datasets)) {
  ds_dir <- file.path(results_root, ds)
  if (!dir.exists(ds_dir)) next
  files <- list.files(ds_dir, pattern = "\\.txt$", full.names = TRUE)
  if (!length(files)) next
  parsed <- lapply(files, extract_nll_from_file)
  for (p in parsed) {
    records[[length(records) + 1L]] <- data.frame(
      dataset = ds,
      N = as.integer(p$N),
      seed = as.integer(p$seed),
      nll = as.numeric(p$nll),
      stringsAsFactors = FALSE
    )
  }
}

if (!length(records)) stop("No result files found across datasets.")

df <- do.call(rbind, records)
df <- df[!is.na(df$N) & !is.na(df$nll) & is.finite(df$nll), , drop = FALSE]
if (!nrow(df)) stop("No usable data extracted.")

if (!any(df$N %in% Ns)) Ns <- sort(unique(df$N))

agg_fun <- function(x) c(mean = mean(x), se = stderr(x), n = length(x))
agg <- aggregate(df$nll, by = list(dataset = df$dataset, N = df$N), FUN = agg_fun)
agg_res <- agg$x
if (is.list(agg_res)) {
  stat_mat <- do.call(rbind, agg_res)
} else {
  stat_mat <- agg_res
}
stats <- data.frame(dataset = agg$dataset,
                    N = agg$N,
                    mean = stat_mat[, "mean"],
                    se = stat_mat[, "se"],
                    n = stat_mat[, "n"],
                    stringsAsFactors = FALSE)
stats$se[is.na(stats$se)] <- 0
stats <- stats[order(stats$dataset, stats$N), ]
rownames(stats) <- NULL

Ns_present <- sort(unique(stats$N))
if (!length(Ns_present)) stop("No N values present in results.")

x_ticks <- Ns[Ns %in% Ns_present]
if (!length(x_ticks)) x_ticks <- Ns_present

se_vals <- ifelse(is.na(stats$se), 0, stats$se)
y_min <- min(stats$mean - se_vals, na.rm = TRUE)
y_max <- max(stats$mean + se_vals, na.rm = TRUE)
if (!is.finite(y_min) || !is.finite(y_max)) {
  y_min <- min(stats$mean, na.rm = TRUE)
  y_max <- max(stats$mean, na.rm = TRUE)
}
if (!is.finite(y_min) || !is.finite(y_max)) stop("Unable to determine y-axis limits")

out_png <- file.path(results_root, "N_sensitivity_all.png")

png(out_png, width = 10, height = 6.7, units = "in", res = 400)
par(mar = c(6, 6, 1, 1))
plot(range(x_ticks), c(y_min, y_max), type = "n",
     xlab = "N", ylab = "NLL (lower is better)", xaxt = "n")
axis(1, at = x_ticks, labels = x_ticks)

colors <- c(gas = "#1b9e77", power = "#d95f02", hepmass = "#7570b3", miniboone = "#e7298a")
pch_seq <- c(gas = 16, power = 17, hepmass = 15, miniboone = 18)
lty_seq <- c(gas = 1, power = 2, hepmass = 3, miniboone = 4)

datasets_present <- intersect(names(datasets), unique(stats$dataset))
for (ds in datasets_present) {
  sm <- stats[stats$dataset == ds, , drop = FALSE]
  sm <- sm[order(sm$N), , drop = FALSE]
  col <- if (!is.null(colors[ds])) colors[ds] else "black"
  pch <- if (!is.null(pch_seq[ds])) pch_seq[ds] else 16
  lty <- if (!is.null(lty_seq[ds])) lty_seq[ds] else 1
  lines(sm$N, sm$mean, type = "b", col = col, pch = pch, lty = lty)
  ok_err <- sm$se > 0 & !is.na(sm$se)
  if (any(ok_err)) {
    arrows(sm$N[ok_err], sm$mean[ok_err] - sm$se[ok_err],
           sm$N[ok_err], sm$mean[ok_err] + sm$se[ok_err],
           angle = 90, code = 3, length = 0.05, col = col)
  }
}

if (length(datasets_present)) {
  legend("topright",
         legend = datasets[datasets_present],
         col = colors[datasets_present],
         pch = pch_seq[datasets_present],
         lty = lty_seq[datasets_present],
         bty = "n", title = "Dataset")
}

dev.off()

message(sprintf("Wrote: %s", out_png))

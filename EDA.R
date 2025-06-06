#' Exploratory Data Analysis for Generated Samples
#'
#' Computes distribution of per-sample parameters and pairwise
#' scatter plots with parameter overlay. Results are written to a PDF.
#'
#' @param X matrix of generated samples
#' @param cfg configuration list as in `gen_samples`
#' @param output_file path to PDF file
#' @return invisibly the path to the created PDF
compute_param_values <- function(X, cfg) {
  stopifnot(is.matrix(X), is.list(cfg))
  N <- nrow(X)
  K <- length(cfg)
  params <- vector("list", K)
  for (k in seq_len(K)) {
    c_k <- cfg[[k]]
    if (is.null(c_k$parm)) next
    prev_names <- paste0("X", seq_len(k - 1))
    # determine names via first sample
    if (k == 1) {
      prev_df <- data.frame()
    } else {
      prev_df <- as.data.frame(as.list(X[1, seq_len(k - 1)]))
      names(prev_df) <- prev_names
    }
    args0 <- c_k$parm(prev_df)
    if (c_k$distr == "gamma" && all(c("shape1", "shape2") %in% names(args0))) {
      args0 <- list(shape = args0$shape1, scale = args0$shape2)
    }
    arg_names <- names(args0)
    mat <- matrix(NA_real_, nrow = N, ncol = length(arg_names))
    colnames(mat) <- arg_names
    for (i in seq_len(N)) {
      if (k == 1) {
        prev_df <- data.frame()
      } else {
        prev_df <- as.data.frame(as.list(X[i, seq_len(k - 1)]))
        names(prev_df) <- prev_names
      }
      args <- c_k$parm(prev_df)
      if (c_k$distr == "gamma" && all(c("shape1", "shape2") %in% names(args))) {
        args <- list(shape = args$shape1, scale = args$shape2)
      }
      args <- lapply(args, function(p) if (!is.finite(p) || p <= 0) 1e-3 else p)
      mat[i, ] <- as.numeric(args[arg_names])
    }
    params[[k]] <- as.data.frame(mat)
  }
  params
}

round_df <- function(df, digits = 3) {
  idx <- vapply(df, is.numeric, logical(1))
  df[idx] <- lapply(df[idx], round, digits = digits)
  df
}

create_EDA_report <- function(X, cfg, output_file = "EDA.pdf",
                              scatter_data = NULL, predict_runtime = NULL,
                              table_normal = NULL, table_perm = NULL, perm = NULL) {
  params <- compute_param_values(X, cfg)
  pdf(output_file)

  if (!is.null(predict_runtime)) {
    plot.new()
    title(main = "Laufzeiten der Vorhersage")
    bullet <- "\u2022"
    lines <- c(
      sprintf("%s Baseline: %.3f s (keine Hyperparameter)", bullet, predict_runtime$baseline),
      sprintf("%s TRTF: %.3f s (ntree=%d, mtry=%d, minsplit=%d, minbucket=%d, maxdepth=%d)",
              bullet, predict_runtime$trtf, 50,
              floor(sqrt(ncol(X) - 1)), 25, 20, 4),
      sprintf("%s KS: %.3f s (bw.nrd0)", bullet, predict_runtime$ks),
      sprintf("%s TTM: %.3f s (lr=1e-2, epochs=200, patience=10)", bullet, predict_runtime$ttm)
    )
    y_pos <- 0.9
    step <- 0.07
    for (i in seq_along(lines)) {
      text(0.05, y_pos - step * (i - 1), lines[i], adj = 0)
    }
  }

  if (!is.null(table_normal)) {
    plot.new()
    title(main = "Normale iteration 1 \u2192 2 \u2192 3 \u2192 4")
    tabn <- round_df(table_normal, digits = 3)
    y_pos <- 0.9
    step <- 0.05
    text(0.05, y_pos, paste(names(tabn), collapse = " | "), adj = 0)
    for (i in seq_len(nrow(tabn))) {
      text(0.05, y_pos - step * i, paste(tabn[i, ], collapse = " | "), adj = 0)
    }
  }

  if (!is.null(table_perm)) {
    plot.new()
    if (!is.null(perm)) {
      perm_title <- paste(perm, collapse = " \u2192 ")
    } else {
      perm_title <- "Permutation"
    }
    title(main = paste0("Permutation ", perm_title))
    tabp <- round_df(table_perm, digits = 3)
    y_pos <- 0.9
    step <- 0.05
    text(0.05, y_pos, paste(names(tabp), collapse = " | "), adj = 0)
    for (i in seq_len(nrow(tabp))) {
      text(0.05, y_pos - step * i, paste(tabp[i, ], collapse = " | "), adj = 0)
    }
  }
  for (k in seq_along(params)) {
    p_df <- params[[k]]
    if (is.null(p_df)) next
    for (nm in names(p_df)) {
      hist(p_df[[nm]], breaks = 30, main = sprintf("%s for X%d", nm, k),
           xlab = nm, col = "grey", border = "white")
    }
  }
  for (k in 2:ncol(X)) {
    p_df <- params[[k]]
    if (is.null(p_df) || ncol(p_df) < 1) next
    color_val <- p_df[[1]]
    param_name <- names(p_df)[1]
    cols <- colorRampPalette(c("blue", "red"))(100)
    idx <- cut(color_val, breaks = 100, labels = FALSE, include.lowest = TRUE)
    plot(X[, k - 1], X[, k], col = cols[idx], pch = 16,
         xlab = paste0("X", k - 1), ylab = paste0("X", k),
         main = sprintf("X%d vs X%d by %s", k, k - 1, param_name))
    legend("topright", legend = c(sprintf("low %0.2f", min(color_val)),
                                 sprintf("high %0.2f", max(color_val))),
           fill = cols[c(1, 100)], bty = "n", title = param_name)
  }

 if (!is.null(scatter_data)) {
    with(scatter_data, {
      plot(ld_base, ld_base, pch = 16, col = "black",
           xlab = "True log-Dichte",
           ylab = "Geschaetzte log-Dichte",
           main = "Originalreihenfolge")
      points(ld_base, ld_trtf, col = "blue", pch = 1)
      points(ld_base, ld_ks,   col = "red",  pch = 2)
      points(ld_base, ld_ttm,  col = "darkgreen", pch = 3)
      abline(a = 0, b = 1)
      legend("topleft", legend = c("true_baseline", "trtf", "ks", "ttm"),
             col = c("black", "blue", "red", "darkgreen"), pch = c(16, 1, 2, 3))

      plot(ld_base_p, ld_base_p, pch = 16, col = "black",
           xlab = "True log-Dichte",
           ylab = "Geschaetzte log-Dichte",
           main = "Permutation")
      points(ld_base_p, ld_trtf_p, col = "blue", pch = 1)
      points(ld_base_p, ld_ks_p,   col = "red",  pch = 2)
      points(ld_base_p, ld_ttm_p,  col = "darkgreen", pch = 3)
      abline(a = 0, b = 1)
      legend("topleft", legend = c("true_baseline", "trtf", "ks", "ttm"),
             col = c("black", "blue", "red", "darkgreen"), pch = c(16, 1, 2, 3))
    })
  }
  
  dev.off()
  invisible(output_file)
}

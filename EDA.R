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

create_EDA_report <- function(X, cfg, output_file = "eda_report.pdf") {
  params <- compute_param_values(X, cfg)
  pdf(output_file)
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
  dev.off()
  invisible(output_file)
}

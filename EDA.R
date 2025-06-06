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

create_EDA_report <- function(X, cfg, output_file = "eda_report.pdf",
                              scatter_data = NULL) {
  pdf(output_file, width = 8.27, height = 11.69)

  if (!is.null(scatter_data)) {
    par(mfrow = c(3, 2))
    with(scatter_data, {
      plot(ld_base, ld_trtf, pch = 16, col = "blue",
           xlab = "True log-Dichte",
           ylab = "TRTF log-Dichte",
           main = "TRTF vs true")
      abline(a = 0, b = 1)

      plot(ld_base_p, ld_trtf_p, pch = 16, col = "blue",
           xlab = "True log-Dichte",
           ylab = "TRTF log-Dichte",
           main = "TRTF vs true (perm)")
      abline(a = 0, b = 1)

      plot(ld_base, ld_ks, pch = 16, col = "red",
           xlab = "True log-Dichte",
           ylab = "KS log-Dichte",
           main = "KS vs true")
      abline(a = 0, b = 1)

      plot(ld_base_p, ld_ks_p, pch = 16, col = "red",
           xlab = "True log-Dichte",
           ylab = "KS log-Dichte",
           main = "KS vs true (perm)")
      abline(a = 0, b = 1)

      plot(ld_base, ld_ttm, pch = 16, col = "darkgreen",
           xlab = "True log-Dichte",
           ylab = "TTM log-Dichte",
           main = "TTM vs true")
      abline(a = 0, b = 1)

      plot(ld_base_p, ld_ttm_p, pch = 16, col = "darkgreen",
           xlab = "True log-Dichte",
           ylab = "TTM log-Dichte",
           main = "TTM vs true (perm)")
      abline(a = 0, b = 1)
    })
  }

  dev.off()
  invisible(output_file)
}

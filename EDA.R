#' Exploratory Data Analysis for Generated Samples
#'
#' Creates tables and scatter plots comparing estimated log-densities with
#' the wahren Werten. Results are written to a PDF.
#'
#' @param X matrix of generated samples
#' @param cfg configuration list as in `gen_samples`
#' @param output_file path to PDF file
#' @return invisibly the path to the created PDF

round_df <- function(df, digits = 3) {
  idx <- vapply(df, is.numeric, logical(1))
  df[idx] <- lapply(df[idx], round, digits = digits)
  df
}

create_EDA_report <- function(X, cfg, output_file = "eda_report.pdf",
                              scatter_data = NULL,
                              table_kbl = NULL) {
  if (!requireNamespace("gridExtra", quietly = TRUE))
    install.packages("gridExtra", repos = "https://cloud.r-project.org")



  make_plot <- function(x, y, ttl) {
    ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
      geom_point(color = "steelblue", size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      labs(title = ttl, x = "True log-Dichte", y = "Geschaetzte log-Dichte") +
      theme_minimal()
  }

  plots <- NULL
  if (!is.null(scatter_data)) {
    with(scatter_data, {
      plots <<- list(
        make_plot(ld_base,   ld_trtf,   "TRTF"),
        make_plot(ld_base_p, ld_trtf_p, "TRTF (Perm.)"),
        make_plot(ld_base,   ld_ks,     "KS"),
        make_plot(ld_base_p, ld_ks_p,   "KS (Perm.)"),
        make_plot(ld_base,   ld_ttm,    "TTM"),
        make_plot(ld_base_p, ld_ttm_p,  "TTM (Perm.)")
      )
    })
  }

  pdf(output_file, width = 8, height = 11, title = "Average logL")
  if (!is.null(table_kbl)) {
    tbl <- gridExtra::tableGrob(
      attr(table_kbl, "tab_data"), rows = NULL,
      theme = gridExtra::ttheme_default(base_size = 9)
    )
    gridExtra::grid.arrange(tbl, top = "Average -logLikelihood")
  }
  if (!is.null(plots))
    gridExtra::grid.arrange(grobs = plots, ncol = 2)
  dev.off()
  invisible(output_file)
}

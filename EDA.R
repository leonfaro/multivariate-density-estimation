#' Exploratory Data Analysis for Generated Samples
#'
#' Creates tables and scatter plots comparing estimated log-densities with
#' the wahren Werten. Results are written to a PDF.
#'
#' @param X matrix of generated samples
#' @param cfg configuration list as in `gen_samples`
#' @param output_file path to PDF file
#' @return invisibly the path to the created PDF

create_param_plots <- function(param_list) {
  plots <- list()
  idx <- 1L
  for (k in seq_along(param_list)) {
    df <- param_list[[k]]
    if (k == 1 || is.null(df)) next
    for (nm in names(df)) {
      plt <- ggplot(df, aes_string(x = nm)) +
        geom_histogram(bins = 30, fill = "steelblue", color = "black") +
        labs(title = paste0("X", k, ": ", nm), x = nm, y = "Haeufigkeit") +
        theme_minimal()
      plots[[idx]] <- plt
      idx <- idx + 1L
    }
  }
  plots
}

create_EDA_report <- function(X, cfg, output_file = "eda_report.pdf",
                              scatter_data = NULL,
                              table_kbl = NULL,
                              param_list = NULL) {
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
  param_plots <- NULL
  if (!is.null(scatter_data)) {
    with(scatter_data, {
      plots <<- list(
        make_plot(ld_base,   ld_trtf,   "TRTF"),
        make_plot(ld_base_p, ld_trtf_p, "TRTF (Perm.)"),
        make_plot(ld_base,   ld_ks,     "KS"),
        make_plot(ld_base_p, ld_ks_p,   "KS (Perm.)")
      )
    })
  }
  if (!is.null(param_list))
    param_plots <- create_param_plots(param_list)

  pdf(output_file, width = 8, height = 11, title = "Average logL")
  if (!is.null(table_kbl)) {
    tab_data <- attr(table_kbl, "tab_data")
    num_cols <- vapply(tab_data, is.numeric, logical(1))
    for (col in names(tab_data)[num_cols]) {
      last_val <- as.numeric(tab_data[nrow(tab_data), col])
      tab_data[-nrow(tab_data), col] <- sprintf("%.3f", as.numeric(tab_data[-nrow(tab_data), col]))
      tab_data[nrow(tab_data), col]  <- sprintf("%.0f", last_val)
    }
    print(tab_data)
  }
  if (!is.null(plots))
    gridExtra::grid.arrange(grobs = plots, ncol = 2)
  if (!is.null(param_plots)) {
    for (pg in param_plots)
      print(pg)
  }
  dev.off()
  invisible(output_file)
}

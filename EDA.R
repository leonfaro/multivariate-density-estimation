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
                              runtime_list = NULL,
                              hyperparam_list = NULL,
                              tab_normal = NULL,
                              tab_perm = NULL,
                              perm_vec = NULL) {
  if (!requireNamespace("knitr", quietly = TRUE))
    install.packages("knitr", repos = "https://cloud.r-project.org")

  pdf(output_file, paper = "a4")

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  mat <- matrix(c(1,1,1,
                  2,2,2,
                  3,3,3,
                  4,5,6,
                  7,8,9),
                byrow = TRUE, ncol = 3)
  layout(mat, heights = c(1, 1, 1, 1.5, 1.5))
  par(mar = c(2, 2, 3, 1))

  if (!is.null(runtime_list)) {
    plot.new()
    title(main = "Vorhersage-Laufzeiten")
    y_pos <- 0.9
    step <- 0.1
    r_names <- names(runtime_list)
    for (i in seq_along(runtime_list)) {
      nm <- r_names[i]
      hp <- if (!is.null(hyperparam_list) && !is.null(hyperparam_list[[nm]]))
        hyperparam_list[[nm]] else ""
      txt <- sprintf("\u2022 %s (%s): %0.2fs", nm, hp, runtime_list[[i]])
      text(0.05, y_pos - step * (i - 1), txt, adj = 0)
    }
    text(0.05, y_pos - step * (length(runtime_list) + 1),
         paste("N =", nrow(X)), adj = 0)
    cfg_text <- paste(sapply(cfg, `[[`, "distr"), collapse = " \u2192 ")
    text(0.05, y_pos - step * (length(runtime_list) + 2),
         paste("config:", cfg_text), adj = 0)
  } else {
    plot.new()
  }

  if (!is.null(tab_normal)) {
    plot.new()
    title(main = "Normale iteration 1 \u2192 2 \u2192 3 \u2192 4")
    tabn <- round_df(tab_normal, digits = 3)
    tbl_lines <- strsplit(knitr::kable(tabn, format = "simple"), "\n")[[1]]
    y_pos <- 0.9
    step <- 0.05
    for (i in seq_along(tbl_lines)) {
      text(0.05, y_pos - step * (i - 1), tbl_lines[i], adj = 0, family = "mono")
    }
  } else {
    plot.new()
  }

  if (!is.null(tab_perm)) {
    plot.new()
    perm_text <- if (is.null(perm_vec)) "" else paste(perm_vec, collapse = " \u2192 ")
    title(main = sprintf("Permutation %s", perm_text))
    tabp <- round_df(tab_perm, digits = 3)
    tbl_lines <- strsplit(knitr::kable(tabp, format = "simple"), "\n")[[1]]
    y_pos <- 0.9
    step <- 0.05
    for (i in seq_along(tbl_lines)) {
      text(0.05, y_pos - step * (i - 1), tbl_lines[i], adj = 0, family = "mono")
    }
  } else {
    plot.new()
  }
  if (!is.null(scatter_data)) {
    with(scatter_data, {
      plot(ld_base, ld_trtf, pch = 16, col = "blue",
           xlab = "True log-Dichte",
           ylab = "Geschaetzte log-Dichte",
           main = "TRTF")
      abline(a = 0, b = 1)

      plot(ld_base_p, ld_trtf_p, pch = 16, col = "blue",
           xlab = "True log-Dichte",
           ylab = "Geschaetzte log-Dichte",
           main = "TRTF (Perm.)")
      abline(a = 0, b = 1)

      plot(ld_base, ld_ks, pch = 16, col = "red",
           xlab = "True log-Dichte",
           ylab = "Geschaetzte log-Dichte",
           main = "KS")
      abline(a = 0, b = 1)

      plot(ld_base_p, ld_ks_p, pch = 16, col = "red",
           xlab = "True log-Dichte",
           ylab = "Geschaetzte log-Dichte",
           main = "KS (Perm.)")
      abline(a = 0, b = 1)

      plot(ld_base, ld_ttm, pch = 16, col = "darkgreen",
           xlab = "True log-Dichte",
           ylab = "Geschaetzte log-Dichte",
           main = "TTM")
      abline(a = 0, b = 1)

      plot(ld_base_p, ld_ttm_p, pch = 16, col = "darkgreen",
           xlab = "True log-Dichte",
           ylab = "Geschaetzte log-Dichte",
           main = "TTM (Perm.)")
      abline(a = 0, b = 1)
    })
    par(mfrow = c(1, 1))
  }
  
  dev.off()
  invisible(output_file)
}

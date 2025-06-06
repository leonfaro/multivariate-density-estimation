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
  pdf(output_file, paper = "a4")

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
  }

  if (!is.null(tab_normal)) {
    plot.new()
    title(main = "Normale iteration 1 \u2192 2 \u2192 3 \u2192 4")
    tabn <- round_df(tab_normal, digits = 3)
    y_pos <- 0.9
    step <- 0.05
    text(0.05, y_pos, paste(names(tabn), collapse = " | "), adj = 0)
    for (i in seq_len(nrow(tabn))) {
      text(0.05, y_pos - step * i, paste(tabn[i, ], collapse = " | "), adj = 0)
    }
  }

  if (!is.null(tab_perm)) {
    plot.new()
    perm_text <- if (is.null(perm_vec)) "" else paste(perm_vec, collapse = " \u2192 ")
    title(main = sprintf("Permutation %s", perm_text))
    tabp <- round_df(tab_perm, digits = 3)
    y_pos <- 0.9
    step <- 0.05
    text(0.05, y_pos, paste(names(tabp), collapse = " | "), adj = 0)
    for (i in seq_len(nrow(tabp))) {
      text(0.05, y_pos - step * i, paste(tabp[i, ], collapse = " | "), adj = 0)
    }
  }



  if (!is.null(scatter_data)) {
    par(mfrow = c(3, 2))
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

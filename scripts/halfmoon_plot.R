#' Fit half-moon models and create contour plots
#'
#' Provides two helper functions for working with the two-moons data
#' splits.

fit_halfmoon_models <- function(S, seed = NULL) {
  seed <- if (!is.null(seed)) as.integer(seed) else if (!is.null(S$meta$seed)) S$meta$seed else 42L
  set.seed(seed)
  cfg <- list(list(distr = "norm"), list(distr = "norm"))
  list(
    true = fit_TRUE(S, cfg),
    trtf = fit_TRTF(S, cfg, seed = seed),
    ttm = trainMarginalMap(S, seed = seed)$S,
    ttm_sep = trainSeparableMap(S, seed = seed)$S,
    ttm_cross = trainCrossTermMap(S, seed = seed)$S
  )
}

compute_limits <- function(S, pad = 0.05) {
  Xall <- rbind(S$X_tr, S$X_val, S$X_te)
  xr <- range(Xall[, 1])
  yr <- range(Xall[, 2])
  dx <- pad * diff(xr)
  dy <- pad * diff(yr)
  list(xlim = xr + c(-dx, dx), ylim = yr + c(-dy, dy))
}

draw_points <- function(S, style = list()) {
  n_all <- nrow(S$X_tr) + nrow(S$X_val) + nrow(S$X_te)
  cex <- if (!is.null(style$cex)) style$cex else min(1.2, sqrt(800 / n_all))
  cols <- if (!is.null(style$cols)) style$cols else c(
    train = rgb(0.2, 0.4, 1, 0.65),
    val = rgb(1, 0.6, 0, 0.5),
    test = rgb(1, 0, 0, 0.7)
  )
  points(S$X_tr, pch = 16, cex = cex, col = cols["train"])
  points(S$X_val, pch = 16, cex = cex, col = cols["val"])
  points(S$X_te, pch = 16, cex = cex * 1.1, col = cols["test"])
  invisible(n_all)
}

clip01 <- function(x) {
  q <- stats::quantile(x, c(0.01, 0.99))
  pmin(pmax(x, q[1]), q[2])
}

.draw_panels <- function(mods, S, grid_n, levels_policy = c("global_quantiles", "per_model")) {
  levels_policy <- match.arg(levels_policy)
  lim <- compute_limits(S, pad = 0.05)
  xseq <- seq(lim$xlim[1], lim$xlim[2], length.out = grid_n)
  yseq <- seq(lim$ylim[1], lim$ylim[2], length.out = grid_n)
  G <- as.matrix(expand.grid(xseq, yseq))
  get_LDj_true <- function(mod_true) {
    cbind(
      .log_density_vec(G[, 1], "norm", mod_true$theta[[1]]),
      .log_density_vec(G[, 2], "norm", mod_true$theta[[2]])
    ) |> rowSums()
  }
  get_LDj <- function(name) {
    if (name == "true") {
      get_LDj_true(mods$true)
    } else {
      as.numeric(predict(mods[[name]], G, "logdensity"))
    }
  }
  panels <- c("true", "trtf", "ttm", "ttm_sep", "ttm_cross")
  LD_list <- setNames(lapply(panels, get_LDj), panels)
  all_finite <- all(sapply(LD_list, function(z) all(is.finite(z))))
  probs <- c(0.9, 0.7, 0.5)
  if (levels_policy == "global_quantiles") {
    all_vals <- unlist(lapply(LD_list, function(z) clip01(z[is.finite(z)])))
    lev <- stats::quantile(all_vals, probs)
    message("Contour levels (global): ", paste(sprintf("%.3f", lev), collapse = ", "))
    lev_list <- rep(list(lev), length(panels))
    names(lev_list) <- panels
    lev_ret <- lev
  } else {
    lev_list <- lapply(LD_list, function(z) {
      stats::quantile(clip01(z[is.finite(z)]), probs)
    })
    message("Contour levels (per_model): ",
            paste(names(lev_list),
                  sapply(lev_list, function(v) paste(sprintf("%.3f", v), collapse = ", ")), collapse = " | "))
    lev_ret <- lev_list
  }
  n_all <- nrow(S$X_tr) + nrow(S$X_val) + nrow(S$X_te)
  style <- list(
    cex = min(1.2, sqrt(800 / n_all)),
    cols = c(train = rgb(0.2, 0.4, 1, 0.65),
             val = rgb(1, 0.6, 0, 0.5),
             test = rgb(1, 0, 0, 0.7))
  )
  op <- par(mfrow = c(2, 3), mar = c(3, 3, 2, 1))
  plot(NA, xlim = lim$xlim, ylim = lim$ylim, xlab = "x1", ylab = "x2", main = "Data")
  grid(col = "gray90", lty = 3)
  draw_points(S, style)
  legend("bottomright", legend = c("Train", "Val", "Test"), pch = 16,
         col = style$cols, pt.cex = c(style$cex, style$cex, style$cex * 1.1), bty = "n")
  for (nm in panels) {
    Z <- matrix(LD_list[[nm]], nrow = length(xseq), ncol = length(yseq))
    plot(NA, xlim = lim$xlim, ylim = lim$ylim, xlab = "x1", ylab = "x2", main = nm)
    grid(col = "gray90", lty = 3)
    contour(xseq, yseq, Z, levels = lev_list[[nm]], add = TRUE, drawlabels = FALSE)
    draw_points(S, style)
    legend("bottomright", legend = c("Train", "Val", "Test"), pch = 16,
           col = style$cols, pt.cex = c(style$cex, style$cex, style$cex * 1.1), bty = "n")
  }
  par(op)
  invisible(list(levels = lev_ret, grid_points = nrow(G), all_finite = all_finite))
}

plot_halfmoon_models <- function(mods, S, grid_n = 200,
                                 levels_policy = c("global_quantiles", "per_model"),
                                 save_png = TRUE, show_plot = TRUE,
                                 seed = NULL) {
  levels_policy <- match.arg(levels_policy)
  if (!is.null(seed)) {
    set.seed(seed)
  } else if (!is.null(S$meta$seed)) {
    set.seed(S$meta$seed)
  }
  if (save_png) {
    dir.create("results", showWarnings = FALSE)
    seed0 <- if (!is.null(S$meta$seed)) S$meta$seed else 0
    f <- sprintf("results/halfmoon_panels_seed%03d.png", seed0)
    png(f, width = 1200, height = 800)
    .draw_panels(mods, S, grid_n, levels_policy)
    dev.off()
    message("Saved: ", f)
  }
  if (show_plot && interactive()) {
    .draw_panels(mods, S, grid_n, levels_policy)
  }
  invisible(TRUE)
}

#' ggplot2 panels for half-moon models
#'
#' Creates a 2x3 panel figure comparing data and model densities using
#' `ggplot2` with filled contours at global quantile levels.
plot_halfmoon_models_gg <- function(mods, S, grid_n = 200,
                                    levels_policy = c("global_quantiles"),
                                    save_png = TRUE, show_plot = FALSE) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  levels_policy <- match.arg(levels_policy)
  lim <- compute_limits(S, pad = 0.05)
  xseq <- seq(lim$xlim[1], lim$xlim[2], length.out = grid_n)
  yseq <- seq(lim$ylim[1], lim$ylim[2], length.out = grid_n)
  G <- expand.grid(x1 = xseq, x2 = yseq)
  get_LDj_true <- function(mod_true) {
    cbind(
      .log_density_vec(G$x1, "norm", mod_true$theta[[1]]),
      .log_density_vec(G$x2, "norm", mod_true$theta[[2]])
    ) |> rowSums()
  }
  get_LDj <- function(name) {
    if (name == "true") get_LDj_true(mods$true)
    else as.numeric(predict(mods[[name]], as.matrix(G), "logdensity"))
  }
  panels <- c("true", "trtf", "ttm", "ttm_sep", "ttm_cross")
  LD_list <- setNames(lapply(panels, get_LDj), panels)
  probs <- c(0.9, 0.7, 0.5)
  all_vals <- unlist(lapply(LD_list, function(z) clip01(z[is.finite(z)])))
  lev <- stats::quantile(all_vals, probs)
  message("Contour levels (global): ", paste(sprintf("%.3f", lev), collapse = ", "))
  df_ld <- do.call(rbind, lapply(panels, function(p) {
    data.frame(G, LD = LD_list[[p]], panel = p)
  }))
  panel_levels <- c("data", panels)
  make_df <- function(X, y, split) {
    data.frame(x1 = X[, 1], x2 = X[, 2], y = y, split = split)
  }
  df_pts <- do.call(rbind, lapply(panel_levels, function(pl) {
    rbind(
      cbind(make_df(S$X_tr, S$y_tr, "Train"), panel = pl),
      cbind(make_df(S$X_val, S$y_val, "Val"), panel = pl),
      cbind(make_df(S$X_te, S$y_te, "Test"), panel = pl)
    )
  }))
  df_ld$panel <- factor(df_ld$panel, levels = panel_levels[-1])
  df_pts$panel <- factor(df_pts$panel, levels = panel_levels)
  p <- ggplot2::ggplot() +
    ggplot2::stat_contour_filled(data = df_ld,
      ggplot2::aes(x1, x2, z = LD, fill = after_stat(level)),
      breaks = lev) +
    ggplot2::geom_point(data = df_pts,
      ggplot2::aes(x1, x2, color = factor(y), shape = split),
      alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("1" = "blue", "2" = "red"), name = "Mond") +
    ggplot2::scale_shape_manual(values = c(Train = 16, Val = 17, Test = 15), name = "Split") +
    ggplot2::coord_fixed(xlim = lim$xlim, ylim = lim$ylim) +
    ggplot2::facet_wrap(~panel, nrow = 2) +
    ggplot2::theme_minimal() +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(title = "log-density"),
      color = ggplot2::guide_legend(title = "Mond"),
      shape = ggplot2::guide_legend(title = "Split")
    )
  if (save_png) {
    dir.create("results", showWarnings = FALSE)
    seed0 <- if (!is.null(S$meta$seed)) S$meta$seed else 0
    f <- sprintf("results/halfmoon_panels_seed%03d.png", seed0)
    ggplot2::ggsave(f, p, width = 12, height = 8, dpi = 100)
    message("Saved: ", f)
  }
  if (show_plot && interactive()) {
    print(p)
  }
  invisible(TRUE)
}


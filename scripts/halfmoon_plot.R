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

plot_halfmoon_models_gg <- function(mods, S, grid_n = 200,
                                    levels_policy = c("global_quantiles"),
                                    save_png = TRUE, show_plot = TRUE) {
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
  all_vals <- unlist(LD_list)
  lev <- stats::quantile(all_vals[is.finite(all_vals)], c(0.9, 0.7, 0.5))
  lev <- sort(lev)
  message("Contour levels (global): ", paste(sprintf("%.3f", lev), collapse = ", "))
  breaks <- c(-Inf, lev, Inf)
  grid_df <- do.call(rbind, lapply(panels, function(nm) {
    data.frame(x1 = G[, 1], x2 = G[, 2], ld = LD_list[[nm]], panel = nm)
  }))
  points_df <- rbind(
    data.frame(x1 = S$X_tr[, 1], x2 = S$X_tr[, 2], y = S$y_tr, split = "train"),
    data.frame(x1 = S$X_val[, 1], x2 = S$X_val[, 2], y = S$y_val, split = "val"),
    data.frame(x1 = S$X_te[, 1], x2 = S$X_te[, 2], y = S$y_te, split = "test")
  )
  panels_all <- c("Data", panels)
  point_panel_df <- do.call(rbind, lapply(panels_all, function(p)
    cbind(points_df, panel = p)))
  grid_df$panel <- factor(grid_df$panel, levels = panels_all)
  point_panel_df$panel <- factor(point_panel_df$panel, levels = panels_all)
  library(ggplot2)
  p <- ggplot() +
    stat_contour_filled(data = grid_df,
                        aes(x1, x2, z = ld, fill = after_stat(level)),
                        breaks = breaks) +
    geom_point(data = point_panel_df,
               aes(x1, x2, color = factor(y), shape = split),
               alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(`1` = "blue", `2` = "red"), name = "Mond") +
    scale_shape_manual(values = c(train = 16, val = 17, test = 15), name = "Split") +
    scale_fill_viridis_d(option = "plasma", name = "Logdichte", drop = FALSE) +
    facet_wrap(~panel, nrow = 2) +
    coord_fixed(xlim = lim$xlim, ylim = lim$ylim) +
    theme_minimal()
  if (save_png) {
    dir.create("results", showWarnings = FALSE)
    seed0 <- if (!is.null(S$meta$seed)) S$meta$seed else 0
    f <- sprintf("results/halfmoon_panels_seed%03d.png", seed0)
    ggsave(f, plot = p, width = 12, height = 8)
    message("Saved: ", f)
  }
  if (show_plot && interactive()) {
    print(p)
  }
  invisible(TRUE)
}


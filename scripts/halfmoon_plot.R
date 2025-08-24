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
    ttm = trainMarginalMap(S)$S,
    ttm_sep = trainSeparableMap(S)$S,
    ttm_cross = trainCrossTermMap(S)$S
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

.draw_panels <- function(mods, S, grid_n, levels_policy) {
  lim <- compute_limits(S, pad = 0.05)
  xseq <- seq(lim$xlim[1], lim$xlim[2], length.out = grid_n)
  yseq <- seq(lim$ylim[1], lim$ylim[2], length.out = grid_n)
  G <- as.matrix(expand.grid(xseq, yseq))
  get_LDj_true <- function(mod_true, G) {
    cbind(
      .log_density_vec(G[, 1], "norm", mod_true$theta[[1]]),
      .log_density_vec(G[, 2], "norm", mod_true$theta[[2]])
    ) |> rowSums()
  }
  eval_panel <- function(name) {
    if (name == "true") {
      LDj <- get_LDj_true(mods$true, G)
    } else {
      LDj <- as.numeric(predict(mods[[name]], G, "logdensity"))
    }
    list(LDj = LDj,
         lev = as.numeric(stats::quantile(LDj, probs = levels_policy, na.rm = TRUE)))
  }
  n_all <- nrow(S$X_tr) + nrow(S$X_val) + nrow(S$X_te)
  style <- list(
    cex = min(1.2, sqrt(800 / n_all)),
    cols = c(train = rgb(0.2, 0.4, 1, 0.65),
             val = rgb(1, 0.6, 0, 0.5),
             test = rgb(1, 0, 0, 0.7))
  )
  panels <- c("true", "trtf", "ttm", "ttm_sep", "ttm_cross")
  op <- par(mfrow = c(2, 3), mar = c(3, 3, 2, 1))
  plot(NA, xlim = lim$xlim, ylim = lim$ylim, xlab = "x1", ylab = "x2", main = "Data")
  grid(col = "gray90", lty = 3)
  draw_points(S, style)
  legend("bottomright", legend = c("Train", "Val", "Test"), pch = 16,
         col = style$cols, pt.cex = c(style$cex, style$cex, style$cex * 1.1), bty = "n")
  for (nm in panels) {
    res <- eval_panel(nm)
    Z <- matrix(res$LDj, nrow = length(xseq), ncol = length(yseq))
    plot(NA, xlim = lim$xlim, ylim = lim$ylim, xlab = "x1", ylab = "x2", main = nm)
    grid(col = "gray90", lty = 3)
    contour(xseq, yseq, Z, levels = res$lev, add = TRUE, drawlabels = FALSE)
    draw_points(S, style)
    legend("bottomright", legend = c("Train", "Val", "Test"), pch = 16,
           col = style$cols, pt.cex = c(style$cex, style$cex, style$cex * 1.1), bty = "n")
  }
  par(op)
  invisible(TRUE)
}

plot_halfmoon_models <- function(mods, S, grid_n = 120,
                                 levels_policy = c(0.9, 0.7, 0.5),
                                 save_png = TRUE, show_plot = TRUE,
                                 seed = NULL) {
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


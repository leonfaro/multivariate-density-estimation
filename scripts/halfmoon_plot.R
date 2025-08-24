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

.draw_panels <- function(mods, S, grid_n, levels_policy) {
  Xall <- rbind(S$X_tr, S$X_val, S$X_te)
  xr <- range(Xall[, 1])
  yr <- range(Xall[, 2])
  dx <- 0.1 * diff(xr)
  dy <- 0.1 * diff(yr)
  xseq <- seq(xr[1] - dx, xr[2] + dx, length.out = grid_n)
  yseq <- seq(yr[1] - dy, yr[2] + dy, length.out = grid_n)
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
  panels <- c("true", "trtf", "ttm", "ttm_sep", "ttm_cross")
  op <- par(mfrow = c(2, 3), mar = c(3, 3, 2, 1))
  plot(S$X_te, pch = 16, cex = 0.35, col = gray(0, 0.25),
       xlab = "x1", ylab = "x2", main = "Data")
  abline(h = 0, v = 0, lty = 3, col = "gray")
  for (nm in panels) {
    res <- eval_panel(nm)
    Z <- matrix(res$LDj, nrow = length(xseq), ncol = length(yseq))
    plot(S$X_te, pch = 16, cex = 0.35, col = gray(0, 0.25),
         xlab = "x1", ylab = "x2", main = nm)
    abline(h = 0, v = 0, lty = 3, col = "gray")
    contour(xseq, yseq, Z, levels = res$lev, add = TRUE, drawlabels = FALSE)
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


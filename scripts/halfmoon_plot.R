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

lse <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

iso_mass_levels <- function(joint, xlim, ylim, grid_side, masses) {
  area <- (diff(xlim) / grid_side) * (diff(ylim) / grid_side)
  logp <- joint + log(area)
  p <- exp(logp - lse(logp))
  ord <- order(joint, decreasing = TRUE)
  cs <- cumsum(p[ord])
  sapply(masses, function(m) joint[ord][which(cs >= m)[1]])
}

eval_density_grid <- function(model, G, xlim, ylim, grid_side, seed,
                              type = "logdensity_by_dim", chunk = 2048,
                              cores = min(10L, parallel::detectCores() - 2L),
                              no_cache = FALSE) {
  stopifnot(is.matrix(G))
  model_sig <- digest::digest(utils::capture.output(str(model)))
  cache_key <- digest::digest(list(model_sig, xlim, ylim, grid_side, seed))
  dir.create("results/cache", showWarnings = FALSE, recursive = TRUE)
  cache_path <- file.path("results/cache", sprintf("moon_%s.rds", cache_key))
  if (file.exists(cache_path) && !no_cache) {
    message("Lade Cache: ", cache_path)
    return(readRDS(cache_path))
  }
  N <- nrow(G)
  idx <- split(seq_len(N), ceiling(seq_len(N) / chunk))
  if (getOption("mde.parallel_active", FALSE)) cores <- 1L
  options(mde.parallel_active = TRUE)
  on.exit(options(mde.parallel_active = FALSE), add = TRUE)
  chunk_fun <- function(ii) {
    Xc <- G[ii, , drop = FALSE]
    pr <- try(predict(model, Xc, type = type, cores = 1L), silent = TRUE)
    if (inherits(pr, "try-error")) pr <- predict(model, Xc, type = type)
    stopifnot(is.matrix(pr), nrow(pr) == nrow(Xc))
    pr
  }
  LD_parts <- parallel::mclapply(idx, chunk_fun, mc.cores = cores)
  LD <- do.call(rbind, LD_parts)
  stopifnot(nrow(LD) == N)
  joint <- rowSums(LD)
  if (any(!is.finite(LD))) stop("NA/Inf in log-dichte")
  joint_ref <- try(predict(model, G, type = "logdensity", cores = 1L), silent = TRUE)
  if (inherits(joint_ref, "try-error")) joint_ref <- predict(model, G, "logdensity")
  stopifnot(length(joint_ref) == N, max(abs(joint - joint_ref)) <= 1e-10)
  obj <- list(LD = LD, joint = joint,
              meta = list(cache_key = cache_key, cache_path = cache_path,
                          chunk = chunk, cores = cores, grid_side = grid_side,
                          xlim = xlim, ylim = ylim, seed = seed,
                          model_signature = model_sig))
  saveRDS(obj, cache_path)
  message("Speichere Cache: ", cache_path)
  obj
}

plot_halfmoon_models <- function(mods, S, grid_side,
                                 chunk = 2048,
                                 cores = min(10L, parallel::detectCores() - 2L),
                                 save_png = TRUE, show_plot = TRUE,
                                 no_cache = FALSE) {
  seed <- if (!is.null(S$meta$seed)) S$meta$seed else 0L
  Xtr <- S$X_tr
  mu <- colMeans(Xtr)
  sdv <- apply(Xtr, 2, sd)
  xlim <- mu[1] + c(-3, 3) * sdv[1]
  ylim <- mu[2] + c(-3, 3) * sdv[2]
  xseq <- seq(xlim[1], xlim[2], length.out = grid_side)
  yseq <- seq(ylim[1], ylim[2], length.out = grid_side)
  G <- as.matrix(expand.grid(xseq, yseq))
  panels <- c("true", "trtf", "ttm_cross")
  evals <- setNames(lapply(panels, function(nm) {
    eval_density_grid(mods[[nm]], G, xlim, ylim, grid_side, seed,
                      chunk = chunk, cores = cores, no_cache = no_cache)
  }), panels)
  LD_list <- lapply(evals, `[[`, "joint")
  lev_list <- lapply(LD_list, iso_mass_levels, xlim = xlim, ylim = ylim,
                     grid_side = grid_side, masses = c(0.5, 0.7, 0.9))
  n_all <- nrow(S$X_tr) + nrow(S$X_val) + nrow(S$X_te)
  style <- list(
    cex = min(1.2, sqrt(800 / n_all)),
    cols = c(train = rgb(0.2, 0.4, 1, 0.65),
             val = rgb(1, 0.6, 0, 0.5),
             test = rgb(1, 0, 0, 0.7))
  )
  op <- par(mfrow = c(1, 3), mar = c(3, 3, 2, 1))
  titles <- c(true = "True", trtf = "TRTF", ttm_cross = "Cross-term")
  for (nm in panels) {
    Z <- matrix(LD_list[[nm]], nrow = grid_side, ncol = grid_side)
    plot(NA, xlim = xlim, ylim = ylim, xlab = "x1", ylab = "x2",
         main = titles[[nm]])
    grid(col = "gray90", lty = 3)
    contour(xseq, yseq, Z, levels = lev_list[[nm]], add = TRUE,
            drawlabels = FALSE)
    draw_points(S, style)
  }
  par(op)
  if (save_png) {
    dir.create("results/moon", showWarnings = FALSE, recursive = TRUE)
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    out_dir <- file.path("results/moon", ts)
    dir.create(out_dir)
    png_file <- file.path(out_dir, sprintf("panels_seed%03d.png", seed))
    png(png_file, width = 900, height = 300)
    op2 <- par(mfrow = c(1, 3), mar = c(3, 3, 2, 1))
    for (nm in panels) {
      Z <- matrix(LD_list[[nm]], nrow = grid_side, ncol = grid_side)
      plot(NA, xlim = xlim, ylim = ylim, xlab = "x1", ylab = "x2",
           main = titles[[nm]])
      grid(col = "gray90", lty = 3)
      contour(xseq, yseq, Z, levels = lev_list[[nm]], add = TRUE,
              drawlabels = FALSE)
      draw_points(S, style)
    }
    par(op2)
    dev.off()
    meta <- list(seed = seed, grid_side = grid_side, chunk = chunk,
                 cores = cores, cache_keys = sapply(evals, function(e) e$meta$cache_key))
    saveRDS(meta, file.path(out_dir, "meta.rds"))
    message("Gespeichert: ", png_file)
  }
  if (show_plot && interactive()) {
    invisible() # plot already drawn
  }
  invisible(list(evals = evals, levels = lev_list,
                 png = if (save_png) png_file else NULL))
}


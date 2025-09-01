#' Fit half-moon models and create contour plots
#'
#' Provides two helper functions for working with the two-moons data
#' splits.

fit_halfmoon_models <- function(S, seed = NULL, order_mode = "as-is") {
  seed <- if (!is.null(seed)) as.integer(seed) else if (!is.null(S$meta$seed)) S$meta$seed else 42L
  set.seed(seed)
  cfg <- list(list(distr = "norm"), list(distr = "norm"))
  M_true <- fit_TRUE(S, cfg)
  class(M_true) <- c("true_marginal_model", class(M_true))
  list(
    true = M_true,
    trtf = fit_TRTF(S, cfg, seed = seed),
    ttm = fit_ttm(S, algo = "marginal",  seed = seed)$S,
    ttm_sep = fit_ttm(S, algo = "separable", seed = seed)$S,
    ttm_cross = fit_ttm(S, algo = "crossterm", seed = seed,
                        deg_g = 2L, df_t = 6L, Q = 16L, lambda = 1e-3, maxit = 50L)$S
  )
}

compute_limits <- function(S, pad = 0.05) {
  Xall <- rbind(S$X_tr, S$X_te)
  xr <- range(Xall[, 1])
  yr <- range(Xall[, 2])
  dx <- pad * diff(xr)
  dy <- pad * diff(yr)
  list(xlim = xr + c(-dx, dx), ylim = yr + c(-dy, dy))
}

#' S3-Predict für das TRUE-Marginalmodell, kompatibel zur Predict-API
#'
#' Erwartet ein Objekt mit Feldern `theta` (Liste) und `config` wie aus `fit_TRUE`.
#' Gibt entweder eine N×K-Matrix (logdensity_by_dim) oder die Zeilensumme (logdensity) zurück.
predict.true_marginal_model <- function(object, newdata,
                                        type = c("logdensity", "logdensity_by_dim"),
                                        cores = 1L, ...) {
  type <- match.arg(type)
  stopifnot(is.matrix(newdata))
  cfg <- object$config
  K <- length(cfg)
  N <- nrow(newdata)
  LD <- do.call(cbind, lapply(seq_len(K), function(k) {
    .log_density_vec(newdata[, k], cfg[[k]]$distr, object$theta[[k]])
  }))
  stopifnot(all(is.finite(LD)), nrow(LD) == N, ncol(LD) == K)
  if (type == "logdensity_by_dim") return(LD)
  rowSums(LD)
}

draw_points <- function(S, style = list(), color_by = c("label", "split")) {
  color_by <- match.arg(color_by)
  n_all <- nrow(S$X_tr) + nrow(S$X_te)
  cex <- if (!is.null(style$cex)) style$cex else min(1.2, sqrt(800 / n_all))
  if (color_by == "label") {
    col_tr <- ifelse(S$y_tr == 1L, "blue", "red")
    col_te <- ifelse(S$y_te == 1L, "blue", "red")
    points(S$X_tr, pch = 16, cex = cex, col = adjustcolor(col_tr, alpha.f = 0.65))
    points(S$X_te, pch = 16, cex = cex * 1.1, col = adjustcolor(col_te, alpha.f = 0.8))
  } else {
    cols <- if (!is.null(style$cols)) style$cols else c(
      train = rgb(0.2, 0.4, 1, 0.65),
      test = rgb(1, 0, 0, 0.7)
    )
    points(S$X_tr, pch = 16, cex = cex, col = cols["train"])
    points(S$X_te, pch = 16, cex = cex * 1.1, col = cols["test"])
  }
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

# Gauss–Legendre für [0,1]
gauss_legendre_01 <- function(n) {
  if (n <= 0 || n != as.integer(n)) stop("n must be positive integer")
  if (n == 1) return(list(nodes = 0.5, weights = 1))
  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  for (k in i) {
    J[k, k + 1] <- b[k]
    J[k + 1, k] <- b[k]
  }
  e <- eigen(J, symmetric = TRUE)
  x <- (e$values + 1) / 2
  w <- (2 * (e$vectors[1, ]^2)) / 2
  list(nodes = x, weights = w)
}

# Wahre Half‑Moon Log‑Dichte (Mischung zweier Bögen, gefaltet mit N(0,σ^2 I))
true_halfmoon_logdensity <- function(G, sigma, c = 0.5, Q = 64,
                                     cache_params = list(xlim = NULL, ylim = NULL, grid_side = NULL)) {
  stopifnot(is.matrix(G), ncol(G) == 2)
  dir.create("results/cache", showWarnings = FALSE, recursive = TRUE)
  key_obj <- list(tag = "true_halfmoon", sigma = sigma, c = c, Q = Q,
                  xlim = cache_params$xlim, ylim = cache_params$ylim, grid_side = cache_params$grid_side)
  cache_key <- digest::digest(key_obj)
  cache_path <- file.path("results/cache", sprintf("moon_true_%s.rds", cache_key))
  if (!is.null(cache_params$grid_side) && !is.na(cache_params$grid_side) && file.exists(cache_path)) {
    message("Lade TRUE-Cache: ", cache_path)
    return(readRDS(cache_path))
  }
  gl <- gauss_legendre_01(Q)
  u <- gl$nodes; w <- gl$weights  # sum(w) == 1
  t_nodes <- pi * u
  m1 <- cbind(cos(t_nodes), sin(t_nodes))
  m2 <- cbind(1 - cos(t_nodes), -sin(t_nodes) + c)
  inv2s2 <- 1 / (2 * sigma^2)
  log_norm_const <- -log(2 * pi * sigma^2)
  N <- nrow(G)
  joint <- numeric(N)
  for (i in seq_len(N)) {
    x <- G[i, ]
    d1 <- rowSums((m1 - matrix(x, nrow = length(t_nodes), ncol = 2, byrow = TRUE))^2)
    d2 <- rowSums((m2 - matrix(x, nrow = length(t_nodes), ncol = 2, byrow = TRUE))^2)
    logphi1 <- log_norm_const - inv2s2 * d1
    logphi2 <- log_norm_const - inv2s2 * d2
    a1 <- log(0.5 * w) + logphi1
    a2 <- log(0.5 * w) + logphi2
    all_a <- c(a1, a2)
    m <- max(all_a)
    joint[i] <- m + log(sum(exp(all_a - m)))
  }
  obj <- list(joint = joint,
              meta = list(cache_key = cache_key, cache_path = cache_path,
                          sigma = sigma, Q = Q,
                          xlim = cache_params$xlim, ylim = cache_params$ylim,
                          grid_side = cache_params$grid_side))
  if (!is.null(cache_params$grid_side) && !is.na(cache_params$grid_side)) saveRDS(obj, cache_path)
  obj
}

# HDR-Niveaus (in log-Dichte-Einheiten)
hdr_levels <- function(joint, xlim, ylim, grid_side, masses) {
  area <- (diff(xlim) / grid_side) * (diff(ylim) / grid_side)
  logp <- joint + log(area)
  m <- max(logp)
  logZ <- m + log(sum(exp(logp - m)))
  p <- exp(logp - logZ)
  ord <- order(joint, decreasing = TRUE)
  cs <- cumsum(p[ord])
  sapply(masses, function(mm) joint[ord][which(cs >= mm)[1]])
}

eval_density_grid <- function(model, G, xlim, ylim, grid_side, seed,
                              type = "logdensity_by_dim", chunk = 2048,
                              cores = min(10L, parallel::detectCores() - 2L),
                              no_cache = FALSE, timeout_sec = NA_integer_,
                              abort_file = NULL) {
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
  idx_list <- split(seq_len(N), ceiling(seq_len(N) / chunk))
  n_chunks <- length(idx_list)
  # Chunk-Checkpoint-Verzeichnis
  chunk_dir <- file.path("results/cache", paste0("chunks_", cache_key))
  dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
  chunk_path <- function(i) file.path(chunk_dir, sprintf("chunk_%05d.rds", i))

  # Starte PSOCK-Cluster nur wenn cores_use > 1, sonst sequentiell
  cores_use <- if (is.finite(cores) && cores >= 1L) as.integer(cores) else 1L
  use_cluster <- cores_use > 1L
  if (use_cluster) {
    cl <- parallel::makeCluster(cores_use, type = "PSOCK")
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, {
      options(mde.parallel_active = TRUE)
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                 MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
      NULL
    })
    parallel::clusterEvalQ(cl, {
      source("00_globals.R")
      source("R/ttm_bases.R"); source("R/ttm_core.R")
      source("R/ttm_marginal.R"); source("R/ttm_separable.R"); source("R/ttm_crossterm.R")
      source("models/trtf_model.R"); source("models/true_model.R")
      source("scripts/halfmoon_plot.R")
      NULL
    })
    parallel::clusterExport(cl, varlist = c("G", "type", "model", "idx_list",
                                            "chunk_path", "abort_file", "timeout_sec"),
                            envir = environment())
  }

  # Entscheide, welche Chunks noch fehlen
  missing_ids <- which(!file.exists(vapply(seq_len(n_chunks), chunk_path, character(1))))
  if (length(missing_ids) > 0L) {
    do_chunk_local <- function(i) {
      if (!is.null(abort_file) && nzchar(abort_file) && file.exists(abort_file)) stop("ABORT: abort_file present")
      Xc <- G[idx_list[[i]], , drop = FALSE]
      pr <- try(predict(model, Xc, type = type), silent = TRUE)
      if (inherits(pr, "try-error")) pr <- predict(model, Xc, type = type)
      stopifnot(is.matrix(pr), nrow(pr) == nrow(Xc))
      jointc <- rowSums(pr)
      joint_ref <- try(predict(model, Xc, type = "logdensity"), silent = TRUE)
      if (inherits(joint_ref, "try-error")) joint_ref <- predict(model, Xc, "logdensity")
      stopifnot(length(joint_ref) == nrow(Xc), max(abs(jointc - joint_ref)) <= 1e-10)
      if (any(!is.finite(pr))) stop("NA/Inf in log-dichte (chunk)")
      saveRDS(list(LD = pr, joint = jointc, id = i), chunk_path(i))
      i
    }
    if (use_cluster) {
      par_fun <- function(i) do_chunk_local(i)
      parallel::parLapply(cl, missing_ids, par_fun)
    } else {
      lapply(missing_ids, do_chunk_local)
    }
  }
  # Reassembling
  parts <- lapply(seq_len(n_chunks), function(i) readRDS(chunk_path(i)))
  # Sicherheits-Check
  stopifnot(all(vapply(parts, function(x) x$id, integer(1)) == seq_len(n_chunks)))
  LD <- do.call(rbind, lapply(parts, `[[`, "LD"))
  stopifnot(nrow(LD) == N)
  joint <- rowSums(LD)
  if (any(!is.finite(LD))) stop("NA/Inf in log-dichte")
  joint_ref <- try(predict(model, G, type = "logdensity", cores = 1L), silent = TRUE)
  if (inherits(joint_ref, "try-error")) joint_ref <- predict(model, G, "logdensity")
  stopifnot(length(joint_ref) == N, max(abs(joint - joint_ref)) <= 1e-10)
  obj <- list(LD = LD, joint = joint,
              meta = list(cache_key = cache_key, cache_path = cache_path,
                          chunk = chunk, cores = cores_use, grid_side = grid_side,
                          xlim = xlim, ylim = ylim, seed = seed,
                          model_signature = model_sig, chunk_dir = chunk_dir,
                          timeout_sec = timeout_sec, abort_file = abort_file))
  saveRDS(obj, cache_path)
  message("Speichere Cache: ", cache_path)
  obj
}

plot_halfmoon_models <- function(mods, S, grid_side,
                                 chunk = 2048,
                                 cores = min(10L, parallel::detectCores() - 2L),
                                 save_png = TRUE, show_plot = TRUE,
                                 no_cache = FALSE, timeout_sec = NA_integer_, abort_file = NULL,
                                 order_mode = "as-is") {
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) || cores < 1L) cores <- 1L
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
  # TRUE-Panel: exakte Log-Dichte gem. Generator
  source("scripts/true_halfmoon_density.R")
  true_eval <- true_logdensity(G, S, Q = 32L)
  evals <- list(true = list(LD = true_eval$by_dim, joint = true_eval$joint, meta = list(tag = "true_halfmoon")))
  for (nm in panels[panels != "true"]) {
    evals[[nm]] <- eval_density_grid(mods[[nm]], G, xlim, ylim, grid_side, seed,
                                     chunk = chunk, cores = cores, no_cache = no_cache,
                                     timeout_sec = timeout_sec, abort_file = abort_file)
  }
  LD_list <- lapply(evals, `[[`, "joint")
  lev_list <- lapply(LD_list, hdr_levels, xlim = xlim, ylim = ylim,
                     grid_side = grid_side, masses = c(0.5, 0.7, 0.9))
  n_all <- nrow(S$X_tr) + nrow(S$X_te)
  style <- list(
    cex = min(1.2, sqrt(800 / n_all)),
    cols = c(train = rgb(0.2, 0.4, 1, 0.65),
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
    draw_points(S, style, color_by = "label")
    legend("topright", legend = c("label 1 (upper) = blue", "label 2 (lower) = red"),
           col = c("blue", "red"), pch = 16, bty = "n", cex = 0.8)
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
      draw_points(S, style, color_by = "label")
      legend("topright", legend = c("label 1 (upper) = blue", "label 2 (lower) = red"),
             col = c("blue", "red"), pch = 16, bty = "n", cex = 0.8)
    }
    par(op2)
    dev.off()
    meta <- list(seed = seed, grid_side = grid_side, chunk = chunk,
                 cores = cores, cache_keys = sapply(evals, function(e) e$meta$cache_key),
                 timeout_sec = timeout_sec, abort_file = abort_file)
    saveRDS(meta, file.path(out_dir, "meta.rds"))
    message("Gespeichert: ", png_file)
  }
  if (show_plot && interactive()) {
    invisible() # plot already drawn
  }
  invisible(list(evals = evals, levels = lev_list,
                 png = if (save_png) png_file else NULL))
}

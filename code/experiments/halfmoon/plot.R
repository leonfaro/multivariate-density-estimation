#' Fit half-moon models and create contour plots
#'
#' Provides two helper functions for working with the two-moons data
#' splits.

if (!exists("locate_repo_loader", inherits = TRUE)) {
  locate_repo_loader <- function() {
    detect_script_path <- function() {
      frames <- sys.frames()
      for (i in rev(seq_along(frames))) {
        fi <- frames[[i]]
        if (!is.null(fi$ofile)) {
          path <- tryCatch(normalizePath(fi$ofile, winslash = "/", mustWork = TRUE),
                          error = function(e) NA_character_)
          if (!is.na(path) && nzchar(path)) return(path)
        }
      }
      args <- commandArgs(trailingOnly = FALSE)
      file_arg <- args[grepl("^--file=", args)]
      if (length(file_arg)) {
        cand <- sub("^--file=", "", file_arg[1])
        path <- tryCatch(normalizePath(cand, winslash = "/", mustWork = TRUE),
                         error = function(e) NA_character_)
        if (!is.na(path) && nzchar(path)) return(path)
      }
      NA_character_
    }

    start_dirs <- character()
    script_path <- detect_script_path()
    if (!is.na(script_path) && nzchar(script_path)) {
      start_dirs <- c(start_dirs, dirname(script_path))
    }
    wd <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = FALSE),
                   error = function(e) getwd())
    start_dirs <- unique(c(start_dirs, wd))
    checked <- character()
    for (start in start_dirs) {
      cur <- start
      repeat {
        cur <- tryCatch(normalizePath(cur, winslash = "/", mustWork = FALSE),
                        error = function(e) cur)
        if (!nzchar(cur) || cur %in% checked) break
        checked <- c(checked, cur)
        cand1 <- file.path(cur, "R", "loader.R")
        if (file.exists(cand1)) {
          return(normalizePath(cand1, winslash = "/", mustWork = TRUE))
        }
        cand2 <- file.path(cur, "code", "R", "loader.R")
        if (file.exists(cand2)) {
          return(normalizePath(cand2, winslash = "/", mustWork = TRUE))
        }
        parent <- dirname(cur)
        if (identical(parent, cur)) break
        cur <- parent
      }
    }
    stop("Could not locate loader.R")
  }
}

fit_halfmoon_models <- function(S, seed = NULL, order_mode = "as-is") {
  seed <- if (!is.null(seed)) as.integer(seed) else if (!is.null(S$meta$seed)) S$meta$seed else 42L
  set.seed(seed)
  cfg <- list(list(distr = "norm"), list(distr = "norm"))
  M_true <- fit_TRUE(S, cfg)
  class(M_true) <- c("true_marginal_model", class(M_true))
  # Ensure copula_np fitter is available without side effects
  if (!exists("fit_copula_np", mode = "function")) {
    root <- if (exists("root_path")) root_path else repo_root()
    src <- file.path(root, "R", "models", "copula_np.R")
    if (file.exists(src)) source(src)
  }
  if (!exists("fit_copula_np", mode = "function")) {
    stop("Copula NP fitter `fit_copula_np` must be available to plot Copula NP panel.")
  }

  timing <- list()
  timing$true <- list(train = 0.0, test = NA_real_)
  timing$true_marg <- list(train = 0.0, test = NA_real_)

  trtf_time <- system.time({
    M_trtf <- tryCatch(fit_TRTF(S, cfg, seed = seed), error = function(e) NULL)
  })[["elapsed"]]
  timing$trtf <- list(train = if (!is.null(M_trtf)) trtf_time else NA_real_, test = NA_real_)

  fit_ttm_m <- fit_ttm(S, algo = "marginal", seed = seed)
  timing$ttm <- list(train = fit_ttm_m$time_train, test = fit_ttm_m$time_pred)

  fit_ttm_s <- fit_ttm(S, algo = "separable", seed = seed)
  timing$ttm_sep <- list(train = fit_ttm_s$time_train, test = fit_ttm_s$time_pred)

  cop_time <- system.time({
    M_cop <- fit_copula_np(S, seed = seed)
  })[["elapsed"]]
  timing$copula_np <- list(train = cop_time, test = NA_real_)

  mods <- list(
    true = M_true,
    true_marg = M_true,
    trtf = M_trtf,
    ttm = fit_ttm_m$S,
    ttm_sep = fit_ttm_s$S,
    copula_np = M_cop
  )
  attr(mods, "timing") <- timing
  mods
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
  # Caching disabled: compute directly without reading/writing disk
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
              meta = list(sigma = sigma, Q = Q,
                          xlim = cache_params$xlim, ylim = cache_params$ylim,
                          grid_side = cache_params$grid_side))
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
                              no_cache = TRUE, timeout_sec = NA_integer_,
                              abort_file = NULL) {
  stopifnot(is.matrix(G))
  loader_path <- locate_repo_loader()
  plot_script_path <- tryCatch(
    normalizePath(file.path(dirname(dirname(loader_path)),
                             "experiments", "halfmoon", "plot.R"),
                  winslash = "/", mustWork = TRUE),
    error = function(e) stop("Could not determine path to halfmoon/plot.R", call. = FALSE)
  )
  model_sig <- digest::digest(utils::capture.output(str(model)))
  # Caching disabled: compute directly without reading/writing disk
  N <- nrow(G)
  idx_list <- split(seq_len(N), ceiling(seq_len(N) / chunk))
  n_chunks <- length(idx_list)
  # No on-disk chunk checkpoints
  chunk_path <- function(i) tempfile(sprintf("chunk_%05d_", i))

  # Starte PSOCK-Cluster nur wenn cores_use > 1, sonst sequentiell
  cores_use <- if (is.finite(cores) && cores >= 1L) as.integer(cores) else 1L
  use_cluster <- cores_use > 1L
  if (use_cluster) {
    cl <- parallel::makeCluster(cores_use, type = "PSOCK")
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterExport(cl, varlist = c("loader_path", "plot_script_path"),
                            envir = environment())
    parallel::clusterEvalQ(cl, {
      options(mde.parallel_active = TRUE)
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                 MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
      NULL
    })
    parallel::clusterEvalQ(cl, {
      if (!exists("initialize_repo")) {
        source(loader_path, chdir = FALSE)
      }
      root <- initialize_repo()
      source(plot_script_path, chdir = FALSE)
      NULL
    })
    parallel::clusterExport(cl, varlist = c("G", "type", "model", "idx_list",
                                            "chunk_path", "abort_file", "timeout_sec"),
                            envir = environment())
  }

  # Entscheide, welche Chunks noch fehlen
  missing_ids <- seq_len(n_chunks)
  if (length(missing_ids) > 0L) {
    do_chunk_local <- function(i) {
      if (!is.null(abort_file) && nzchar(abort_file) && file.exists(abort_file)) stop("ABORT: abort_file present")
      Xc <- G[idx_list[[i]], , drop = FALSE]
      pr <- try(predict(model, Xc, type = type), silent = TRUE)
      if (inherits(pr, "try-error")) pr <- predict(model, Xc, type = type)
      stopifnot(is.matrix(pr), nrow(pr) == nrow(Xc))
      jointc <- rowSums(pr)
      if (any(!is.finite(pr))) stop("NA/Inf in log-dichte (chunk)")
      list(LD = pr, joint = jointc, id = i)
    }
    if (use_cluster) {
      parts <- parallel::parLapply(cl, missing_ids, do_chunk_local)
    } else {
      parts <- lapply(missing_ids, do_chunk_local)
    }
  }
  # Reassembling from in-memory parts
  if (!exists("parts", inherits = FALSE) || is.null(parts)) parts <- lapply(missing_ids, do_chunk_local)
  # Sicherheits-Check
  stopifnot(all(vapply(parts, function(x) x$id, integer(1)) == seq_len(n_chunks)))
  LD <- do.call(rbind, lapply(parts, `[[`, "LD"))
  stopifnot(nrow(LD) == N)
  joint <- rowSums(LD)
  if (any(!is.finite(LD))) stop("NA/Inf in log-dichte")
  obj <- list(LD = LD, joint = joint,
              meta = list(chunk = chunk, cores = cores_use, grid_side = grid_side,
                          xlim = xlim, ylim = ylim, seed = seed,
                          model_signature = model_sig,
                          timeout_sec = timeout_sec, abort_file = abort_file))
  obj
}

plot_halfmoon_models <- function(mods, S, grid_side,
                                 chunk = 1024,
                                 cores = min(10L, parallel::detectCores() - 2L),
                                 save_png = TRUE, show_plot = TRUE,
                                 no_cache = FALSE, timeout_sec = NA_integer_, abort_file = NULL,
                                 order_mode = "as-is",
                                 output_dir = file.path(root_path, "experiments", "halfmoon", "results"),
                                 file_stub = NULL) {
  start_time <- proc.time()
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
  panels <- c("true", "true_marg", "trtf", "ttm", "ttm_sep", "copula_np")
  eval_timer_start <- proc.time()
  # TRUE-Panel: exakte Log-Dichte gem. Generator
  source(file.path(root_path, "experiments", "halfmoon", "true_density.R"))
  true_eval <- true_logdensity(G, S, Q = 16L)
  evals <- list(true = list(LD = true_eval$by_dim, joint = true_eval$joint, meta = list(tag = "true_halfmoon")))
  # Helper: unconditional copula_np density over grid by mixing class components
  eval_copula_grid <- function(mod_cnp, G) {
    ld <- predict(mod_cnp, G, type = "logdensity_by_dim")
    lj <- rowSums(ld)
    list(LD = ld, joint = lj, meta = list(tag = "copula_np_grid"))
  }
  missing_models <- setdiff(panels, "true")[vapply(setdiff(panels, "true"), function(nm) is.null(mods[[nm]]), logical(1L))]
  if (length(missing_models) > 0L) {
    stop("Missing fitted models for panels: ", paste(missing_models, collapse = ", "))
  }
  for (nm in panels[panels != "true"]) {
    if (nm == "copula_np") {
      evals[[nm]] <- eval_copula_grid(mods[[nm]], G)
    } else {
      evals[[nm]] <- eval_density_grid(mods[[nm]], G, xlim, ylim, grid_side, seed,
                                       chunk = chunk, cores = cores, no_cache = no_cache,
                                       timeout_sec = timeout_sec, abort_file = abort_file)
    }
  }
  eval_elapsed <- (proc.time() - eval_timer_start)[["elapsed"]]
  LD_list <- lapply(evals[panels], `[[`, "joint")
  lev_list <- lapply(LD_list, hdr_levels, xlim = xlim, ylim = ylim,
                     grid_side = grid_side, masses = c(0.5, 0.7, 0.9))
  n_all <- nrow(S$X_tr) + nrow(S$X_te)
  style <- list(
    cex = min(1.2, sqrt(800 / n_all)),
    cols = c(train = rgb(0.2, 0.4, 1, 0.65),
             test = rgb(1, 0, 0, 0.7))
  )
  n_panels <- length(panels)
  ncol <- min(3L, n_panels)
  nrow <- ceiling(n_panels / ncol)
  op <- par(mfrow = c(nrow, ncol), mar = c(3, 3, 2, 1))
  titles <- c(true = "True (Joint)", true_marg = "True (Marg)", trtf = "TRTF",
              ttm = "TTM Marginal", ttm_sep = "TTM Separable", copula_np = "Copula NP")
  for (nm in panels) {
    Z <- matrix(LD_list[[nm]], nrow = grid_side, ncol = grid_side)
    plot(NA, xlim = xlim, ylim = ylim, xlab = "x1", ylab = "x2",
         main = titles[[nm]])
    grid(col = "gray90", lty = 3)
    contour(xseq, yseq, Z, levels = lev_list[[nm]], add = TRUE,
            drawlabels = FALSE)
    draw_points(S, style, color_by = "label")
  }
  par(op)
  png_elapsed <- NA_real_
  png_file <- NULL
  if (save_png) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    if (is.null(file_stub) || !nzchar(file_stub)) {
      file_stub <- sprintf("halfmoon_panels_seed%03d", seed)
    }
    png_file <- file.path(output_dir, paste0(file_stub, ".png"))
    png_timer_start <- proc.time()
    png(png_file, width = 8, height = 6, units = "in", res = 400)
    op2 <- par(mfrow = c(nrow, ncol), mar = c(3, 3, 2, 1))
    for (nm in panels) {
      Z <- matrix(LD_list[[nm]], nrow = grid_side, ncol = grid_side)
      plot(NA, xlim = xlim, ylim = ylim, xlab = "x1", ylab = "x2",
           main = titles[[nm]])
      grid(col = "gray90", lty = 3)
      contour(xseq, yseq, Z, levels = lev_list[[nm]], add = TRUE,
              drawlabels = FALSE)
      draw_points(S, style, color_by = "label")
    }
    par(op2)
    dev.off()
    png_elapsed <- (proc.time() - png_timer_start)[["elapsed"]]
    meta <- list(seed = seed, grid_side = grid_side, chunk = chunk,
                 cores = cores, output_dir = output_dir,
                 timeout_sec = timeout_sec, abort_file = abort_file)
    saveRDS(meta, file.path(output_dir, paste0(file_stub, "_meta.rds")))
    message("Gespeichert: ", png_file)
  }
  if (show_plot && interactive()) {
    invisible() # plot already drawn
  }
  total_elapsed <- (proc.time() - start_time)[["elapsed"]]
  invisible(list(evals = evals, levels = lev_list,
                 png = if (save_png) png_file else NULL,
                 eval_elapsed = eval_elapsed,
                 plot_elapsed = total_elapsed,
                 png_elapsed = png_elapsed))
}

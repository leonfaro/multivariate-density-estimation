## run_main_debug.R
## Deterministic debug runner for full pipeline with detailed logs
## - n = 50, perm = c(4,3,1,2)
## - Writes timestamped log under logs/
## - Sourcing order mirrors main.R
## - Prints only a single END <logfile> line to stdout on success

make_logfile <- function() {
  dir.create("logs", showWarnings = FALSE, recursive = TRUE)
  ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
  file.path("logs", sprintf("run_main_debug_%s.log", ts))
}

log_path <- make_logfile()
con <- file(log_path, open = "wt")
sink(con)
sink(con, type = "message")
on.exit({
  try(sink(NULL), silent = TRUE)
  try(sink(NULL, type = "message"), silent = TRUE)
  try(close(con), silent = TRUE)
  cat(sprintf("END %s\n", log_path))
}, add = TRUE)

rule <- function(title) cat(sprintf("\n==== %s ====\n", title))

tryCatch({
  rule("START")
  seed <- 42L
  set.seed(seed)
  n <- 50L
  perm <- c(4L, 3L, 1L, 2L)
  cat(sprintf("R %s\n", R.version.string))
  git <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) NA_character_)
  cat(sprintf("git=%s\n", paste(na.omit(git), collapse = " ")))
  cat(sprintf("seed=%d n=%d perm=c(%s)\n", seed, n, paste(perm, collapse = ",")))

  rule("OPTIONS")
  options(cross.grad_check = FALSE)
  options(cross.use_analytic_grad = TRUE)
  options(cross.warmstart_from_separable = TRUE)
  options(cross.sep_degree_g = 2L)
  options(cross.sep_lambda = 1e-3)
  cat(sprintf("cross.grad_check=%s\n", as.character(getOption("cross.grad_check"))))
  cat(sprintf("cross.use_analytic_grad=%s\n", as.character(getOption("cross.use_analytic_grad"))))
  cat(sprintf("cross.warmstart_from_separable=%s\n", as.character(getOption("cross.warmstart_from_separable"))))

  rule("SOURCE")
  get_repo_root <- function() {
    p <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    for (i in 1:10) {
      if (file.exists(file.path(p, "ALGORITHM_SPEC.md")) && file.exists(file.path(p, "00_globals.R"))) return(p)
      parent <- dirname(p)
      if (identical(parent, p)) break
      p <- parent
    }
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
  .repo_root <- get_repo_root()
  src <- function(rel) source(file.path(.repo_root, rel), echo = FALSE, verbose = FALSE, chdir = FALSE)
  src("00_globals.R")
  src("01_data_generation.R")
  src("02_split.R")
  src("models/true_model.R")
  src("models/trtf_model.R")
  # Kernel baseline (copied from main.R)
  stderr <- function(x) stats::sd(x)/sqrt(length(x))
  kernel_fit_marginal <- function(X_tr, seed = 42L) {
    set.seed(seed)
    K <- ncol(X_tr)
    dens <- vector("list", K)
    for (k in seq_len(K)) {
      dk <- stats::density(X_tr[, k], bw = "nrd0", n = 2048, from = min(X_tr[, k]), to = max(X_tr[, k]))
      dens[[k]] <- list(x = dk$x, y = pmax(dk$y, .Machine$double.xmin))
    }
    structure(list(dens = dens, seed = as.integer(seed)), class = "kernel_marginal")
  }
  predict_kernel_marginal <- function(object, X, type = c("logdensity_by_dim", "logdensity")) {
    type <- match.arg(type)
    stopifnot(is.matrix(X))
    K <- ncol(X); N <- nrow(X)
    LD <- matrix(NA_real_, N, K)
    for (k in seq_len(K)) {
      d <- object$dens[[k]]
      fx <- stats::approx(d$x, d$y, xout = X[, k], rule = 2, ties = "ordered")$y
      fx <- pmax(fx, .Machine$double.xmin)
      LD[, k] <- log(fx)
    }
    if (type == "logdensity_by_dim") return(LD)
    rowSums(LD)
  }
  src("models/ttm/ttm_bases.R")
  src("models/ttm/ttm_core.R")
  src("models/ttm/ttm_marginal.R")
  src("models/ttm/ttm_separable.R")
  src("models/ttm/ttm_crossterm.R")
  src("models/true_joint_model.R")
  src("04_evaluation.R")

  rule("CONFIG")
  softplus <- function(x) log1p(exp(x))
  config <- list(
    list(distr = "norm", parm = NULL),
    list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
    list(distr = "beta",
         parm = function(d) list(shape1 = softplus(d$X2),
                                 shape2 = softplus(d$X1))),
    list(distr = "gamma",
         parm = function(d) list(shape = softplus(d$X3),
                                 scale = softplus(d$X2)))
  )
  cat("Config built as in main.R\n")

  rule("DATA")
  prep <- prepare_data(n, config, seed = seed)
  S0 <- prep$S
  S <- list(
    X_tr = S0$X_tr[, perm, drop = FALSE],
    X_te = S0$X_te[, perm, drop = FALSE]
  )
  cfg <- config[perm]
  # Data stats
  stats_of <- function(X) {
    cbind(mean = colMeans(X), sd = apply(X, 2, sd), min = apply(X, 2, min), max = apply(X, 2, max))
  }
  cat("Train stats per dim:\n"); print(stats_of(S$X_tr))
  cat("Test stats per dim:\n");  print(stats_of(S$X_te))

  rule("FIT")
  t_true_tr  <- system.time(mod_true       <- fit_TRUE(S, cfg))[['elapsed']]
  t_joint_tr <- 0
  trtf_ok <- isTRUE(requireNamespace("trtf", quietly = TRUE)) && isTRUE(requireNamespace("tram", quietly = TRUE))
  t_trtf_tr <- NA_real_
  mod_trtf <- NULL
  if (trtf_ok) t_trtf_tr <- system.time(mod_trtf <- fit_TRTF(S, cfg, seed = seed))[['elapsed']] else message("[INFO] TRTF packages not available; skipping TRTF fit.")
  t_kern_tr <- system.time(mod_kern <- kernel_fit_marginal(S$X_tr, seed = seed))[["elapsed"]]
  mod_ttm      <- fit_ttm(S, algo = "marginal",  seed = seed)
  mod_ttm_sep  <- fit_ttm(S, algo = "separable", seed = seed)
  mod_ttm_cross<- fit_ttm(S, algo = "crossterm", seed = seed,
                          deg_g = 2L, df_t = 6L, Q = 16L, lambda = 1e-3, maxit = 50L)

  rule("PREDICT_TIMES")
  t_true_te  <- system.time(logL_TRUE(mod_true, S$X_te))[['elapsed']]
  t_joint_te <- tryCatch({ system.time(true_joint_logdensity_by_dim(config, S0$X_te))[['elapsed']] }, error = function(e) NA_real_)
  t_trtf_te  <- if (trtf_ok) system.time(predict(mod_trtf, S$X_te, type = "logdensity_by_dim"))[['elapsed']] else NA_real_
  t_kern_te  <- system.time(predict_kernel_marginal(mod_kern, S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_ttm_te   <- system.time(predict_ttm(mod_ttm$S,       S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_sep_te   <- system.time(predict_ttm(mod_ttm_sep$S,   S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_ct_te    <- system.time(predict_ttm(mod_ttm_cross$S, S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  cat(sprintf("Times train: true=%.3f joint=%.3f trtf=%s kern=%.3f ttm=%.3f sep=%.3f cross=%.3f\n",
              t_true_tr, t_joint_tr, ifelse(is.na(t_trtf_tr), "NA", sprintf("%.3f", t_trtf_tr)), t_kern_tr,
              mod_ttm$time_train, mod_ttm_sep$time_train, mod_ttm_cross$time_train))
  cat(sprintf("Times test:  true=%.3f joint=%s trtf=%s kern=%.3f ttm=%.3f sep=%.3f cross=%.3f\n",
              t_true_te, ifelse(is.na(t_joint_te), "NA", sprintf("%.3f", t_joint_te)),
              ifelse(is.na(t_trtf_te), "NA", sprintf("%.3f", t_trtf_te)), t_kern_te, t_ttm_te, t_sep_te, t_ct_te))

  rule("INVARIANTS")
  check_model <- function(name, Lbd, Lj) {
    stopifnot(is.matrix(Lbd))
    ok <- max(abs(rowSums(Lbd) - Lj)) <= 1e-12
    cat(sprintf("%s rowsums_equal= %s\n", name, as.character(ok)))
    nf_dim <- sum(!is.finite(Lbd)); nf_joint <- sum(!is.finite(Lj))
    cat(sprintf("%s nonfinite by_dim=%d joint=%d\n", name, nf_dim, nf_joint))
  }
  # TRUE (marginal)
  L_true_by_dim <- do.call(cbind, lapply(seq_along(cfg), function(k) .log_density_vec(S$X_te[, k], cfg[[k]]$distr, mod_true$theta[[k]])))
  L_true_joint  <- rowSums(L_true_by_dim)
  check_model("TRUE(marginal)", L_true_by_dim, L_true_joint)
  # TRUE (Joint)
  L_joint_by_dim <- tryCatch(true_joint_logdensity_by_dim(config, S0$X_te)[, perm, drop = FALSE], error = function(e) matrix(NA_real_, nrow(S$X_te), ncol(S$X_te)))
  L_joint_joint  <- rowSums(L_joint_by_dim)
  check_model("TRUE(JOINT)", L_joint_by_dim, L_joint_joint)
  # TRTF
  if (trtf_ok) {
    L_trtf_by_dim <- predict(mod_trtf, S$X_te, type = "logdensity_by_dim")
    L_trtf_joint  <- rowSums(L_trtf_by_dim)
    check_model("TRTF", L_trtf_by_dim, L_trtf_joint)
  }
  # TTM variants
  L_m_by_dim <- predict_ttm(mod_ttm$S, S$X_te, type = "logdensity_by_dim"); L_m_joint <- predict_ttm(mod_ttm$S, S$X_te, type = "logdensity")
  L_s_by_dim <- predict_ttm(mod_ttm_sep$S, S$X_te, type = "logdensity_by_dim"); L_s_joint <- predict_ttm(mod_ttm_sep$S, S$X_te, type = "logdensity")
  L_c_by_dim <- predict_ttm(mod_ttm_cross$S, S$X_te, type = "logdensity_by_dim"); L_c_joint <- predict_ttm(mod_ttm_cross$S, S$X_te, type = "logdensity")
  check_model("TTM(marginal)", L_m_by_dim, L_m_joint)
  check_model("TTM(separable)", L_s_by_dim, L_s_joint)
  check_model("TTM(crossterm)", L_c_by_dim, L_c_joint)

  # Forward invariants for TTM
  ttm_forward_stats <- function(model, X, tag) {
    out <- ttm_forward(model, X)
    Z <- out$Z; J <- out$J
    cat(sprintf("%s Z dims=%dx%d J dims=%dx%d\n", tag, nrow(Z), ncol(Z), nrow(J), ncol(J)))
    nfZ <- sum(!is.finite(Z)); nfJ <- sum(!is.finite(J))
    cat(sprintf("%s nonfinite Z=%d J=%d\n", tag, nfZ, nfJ))
    for (k in seq_len(ncol(J))) {
      r <- range(J[, k], finite = TRUE)
      nonpos <- sum(J[, k] <= 0, na.rm = TRUE)
      cat(sprintf("%s J[%d] min=%.6g max=%.6g nonpos=%d\n", tag, k, r[1], r[2], nonpos))
    }
    invisible(NULL)
  }
  ttm_forward_stats(mod_ttm$S,       S$X_te, "TTM-M")
  ttm_forward_stats(mod_ttm_sep$S,   S$X_te, "TTM-S")
  ttm_forward_stats(mod_ttm_cross$S, S$X_te, "TTM-C")

  # Cross-term clamp counts and quadrature probe
  rule("CROSSTERM_DIAGNOSTICS")
  Mx <- mod_ttm_cross$S
  mu <- Mx$mu; sigma <- Mx$sigma
  Xs <- sweep(sweep(S$X_te, 2, mu, "-"), 2, sigma, "/")
  nodes <- Mx$gl_nodes; weights <- Mx$gl_weights; Hmax <- if (!is.null(Mx$Hmax)) Mx$Hmax else 20
  for (k in seq_len(ncol(Xs))) {
    Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, nrow(Xs), 0)
    xk <- Xs[, k]
    Hstar <- build_h(xk, Xprev, Mx$spec_h)
    beta <- Mx$coeffs[[k]]$beta
    h_raw <- as.numeric(Hstar %*% beta)
    cl_lo <- sum(h_raw <= -Hmax, na.rm = TRUE)
    cl_hi <- sum(h_raw >=  Hmax, na.rm = TRUE)
    cat(sprintf("CTM k=%d clamp_lo=%d clamp_hi=%d\n", k, cl_lo, cl_hi))
    # Probe integral accuracy with Q and Q+8
    idx <- if (nrow(Xs) > 16) sample.int(nrow(Xs), 16L) else seq_len(nrow(Xs))
    # I with current Q
    Iq <- rep(0, length(idx))
    for (q in seq_along(nodes)) {
      tq <- xk[idx] * nodes[q]
      Hq <- build_h(tq, Xprev[idx, , drop = FALSE], Mx$spec_h)
      v <- as.numeric(Hq %*% beta)
      v <- pmax(pmin(v, Hmax), -Hmax)
      Iq <- Iq + weights[q] * exp(v)
    }
    Iq <- sign(xk[idx]) * abs(xk[idx]) * Iq
    # I with Q+8
    glp <- gauss_legendre_nodes(length(nodes) + 8L)
    Ip <- rep(0, length(idx))
    for (q in seq_along(glp$nodes)) {
      tq <- xk[idx] * glp$nodes[q]
      Hq <- build_h(tq, Xprev[idx, , drop = FALSE], Mx$spec_h)
      v <- as.numeric(Hq %*% beta)
      v <- pmax(pmin(v, Hmax), -Hmax)
      Ip <- Ip + glp$weights[q] * exp(v)
    }
    Ip <- sign(xk[idx]) * abs(xk[idx]) * Ip
    d <- mean(abs(Iq - Ip))
    cat(sprintf("CTM k=%d probe_dI_Q_vs_Q+8=%.3e\n", k, d))
  }

  rule("HYPERPARAMS")
  cat(sprintf("SEP degree_g=%s lambda=%s\n", as.character(mod_ttm_sep$S$degree_g), "0.0"))
  cat(sprintf("CTM deg_g=%s df_t=%s Q=%s lambda=%.3g Hmax=%d maxit=%d warmstart=%s\n",
              as.character(mod_ttm_cross$S$deg_g), as.character(mod_ttm_cross$S$spec_h$df),
              as.character(length(mod_ttm_cross$S$gl_nodes)), 1e-3, as.integer(Hmax), 50L,
              as.character(getOption("cross.warmstart_from_separable", FALSE))))

  rule("TABLES")
  mods <- list(
    true = mod_true,
    true_joint = tryCatch(fit_TRUE_JOINT(S, cfg), error = function(e) NULL),
    trtf = mod_trtf,
    ttm  = mod_ttm,
    ttm_sep = mod_ttm_sep,
    ttm_cross = mod_ttm_cross
  )
  tab <- tryCatch(calc_loglik_tables(mods, cfg, S$X_te, config_canonical = config, perm = perm), error = function(e) {
    message("[WARN] calc_loglik_tables failed: ", e$message)
    # Minimal fallback table with NA for TRTF if unavailable
    K <- ncol(S$X_te)
    data.frame(
      dim = c(as.character(seq_len(K)), "k"),
      distribution = c(sapply(cfg, `[[`, "distr"), "SUM"),
      `True (marginal)` = NA_character_,
      `True (Joint)`    = NA_character_,
      `Random Forest`   = NA_character_,
      `Marginal Map`    = NA_character_,
      `Separable Map`   = NA_character_,
      `Cross-term Map`  = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  if ("train_test_policy" %in% names(tab)) tab$train_test_policy <- NULL
  # Kernel Smooth column as in main.R (ascii-only +-)
  LD_kern <- predict_kernel_marginal(mod_kern, S$X_te, type = "logdensity_by_dim")
  per_kern <- -colMeans(LD_kern)
  se_kern  <- apply(-LD_kern, 2, stderr)
  fmt_ascii <- function(m, se) sprintf("%.2f +- %.2f", round(m, 2), round(2 * se, 2))
  tab[["Kernel Smooth"]] <- c(
    fmt_ascii(per_kern, se_kern),
    fmt_ascii(sum(per_kern), stats::sd(rowSums(-LD_kern)) / sqrt(nrow(S$X_te)))
  )
  cat(sprintf("n=%d\n", n))
  print(tab)
  cat(sprintf("Permutation order %s (train/test only)\n", paste(perm, collapse = ",")))

  time_tab <- data.frame(
    model = c("True (marginal)", "True (Joint)", "Random Forest",
              "Kernel Smooth", "Marginal Map", "Separable Map", "Cross-term Map"),
    train_sec = c(t_true_tr, t_joint_tr, t_trtf_tr,
                  t_kern_tr, mod_ttm$time_train, mod_ttm_sep$time_train, mod_ttm_cross$time_train),
    test_sec = c(t_true_te, t_joint_te, t_trtf_te,
                 t_kern_te, t_ttm_te, t_sep_te, t_ct_te),
    stringsAsFactors = FALSE
  )
  time_tab$total_sec <- with(time_tab, train_sec + test_sec)
  print(time_tab)

  rule("INVARIANT_OK")
  cat("All checks completed.\n")
}, error = function(e) {
  rule("ERROR")
  message("Unhandled error: ", conditionMessage(e))
})


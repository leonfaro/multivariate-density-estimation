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
.repo_root <- if (exists("root_path")) root_path else get_repo_root()
source(file.path(.repo_root, "00_globals.R"))
if (!exists("%||%")) "%||%" <- function(a, b) if (is.null(a)) b else a
source(file.path(.repo_root, "01_data_generation.R"))
source(file.path(.repo_root, "02_split.R"))
source(file.path(.repo_root, "models/true_model.R"))
source(file.path(.repo_root, "models/trtf_model.R"))
# Optional: kernel smoothing marginal baseline (product KDE) helpers
# Deterministic and light-weight; used only for reporting in main().
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
# Use new modular TTM implementation exclusively
source(file.path(.repo_root, "models/ttm/ttm_bases.R"))
source(file.path(.repo_root, "models/ttm/ttm_core.R"))
source(file.path(.repo_root, "models/ttm/ttm_marginal.R"))
source(file.path(.repo_root, "models/ttm/ttm_separable.R"))
source(file.path(.repo_root, "models/ttm/ttm_crossterm.R"))
source(file.path(.repo_root, "models/true_joint_model.R"))
source(file.path(.repo_root, "04_evaluation.R"))

perm <- c(1, 2, 3, 4)
n <- 50
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

#' Starte die komplette Analyse
#' 
#' Es werden zwei Log-Likelihood-Tabellen ausgegeben: einmal für die
#' normale Reihenfolge der Dimensionen und einmal für eine vorgegebene
#' Permutation.
#' 
#' @export
main <- function() {
  seed <- as.integer(Sys.getenv("SEED", 42))
  set.seed(seed)
  prep <- prepare_data(n, config, seed = seed)
  S0 <- prep$S
  S <- list(
    X_tr  = S0$X_tr[, perm, drop = FALSE],
    X_te  = S0$X_te[, perm, drop = FALSE]
  )
  cfg <- config[perm]
  
  t_true_tr  <- system.time(mod_true      <- fit_TRUE(S, cfg))[['elapsed']]
  t_joint_tr <- 0
  mod_true_joint <- tryCatch({
    fit_TRUE_JOINT(S, cfg)
  }, error = function(e) {
    message("[WARN] Skipping True (Joint) for this permutation: ", e$message)
    NULL
  })
  t_trtf_tr  <- system.time(mod_trtf      <- fit_TRTF(S, cfg, seed = seed))[['elapsed']]
  # Kernel smoothing marginal baseline (independent per-dim KDE)
  t_kern_tr  <- system.time(mod_kern      <- kernel_fit_marginal(S$X_tr, seed = seed))[["elapsed"]]
  # New TTM fits (maps-from-samples, forward-KL)
  mod_ttm      <- fit_ttm(S, algo = "marginal",  seed = seed);  t_ttm_tr <- mod_ttm$time_train
  mod_ttm_sep  <- fit_ttm(S, algo = "separable", seed = seed);  t_sep_tr <- mod_ttm_sep$time_train
  # Cross-term with moderate expressiveness and accuracy
  mod_ttm_cross<- fit_ttm(S, algo = "crossterm", seed = seed,
                          deg_g = 2L, df_t = 6L, Q = 16L, lambda = 1e-3, maxit = 50L);
  t_ct_tr  <- mod_ttm_cross$time_train

  t_true_te  <- system.time(logL_TRUE(mod_true, S$X_te))[['elapsed']]
  # Compute True (Joint) timing in canonical order (independent of permutation)
  t_joint_te <- tryCatch({
    system.time(true_joint_logdensity_by_dim(config, S0$X_te))[['elapsed']]
  }, error = function(e) NA_real_)
  t_trtf_te  <- system.time(predict(mod_trtf, S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_kern_te  <- system.time(predict_kernel_marginal(mod_kern, S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_ttm_te   <- system.time(predict_ttm(mod_ttm$S,       S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_sep_te   <- system.time(predict_ttm(mod_ttm_sep$S,   S$X_te, type = "logdensity_by_dim"))[['elapsed']]
  t_ct_te    <- system.time(predict_ttm(mod_ttm_cross$S, S$X_te, type = "logdensity_by_dim"))[['elapsed']]

  mods <- list(
    true = mod_true,
    true_joint = mod_true_joint,
    trtf = mod_trtf,
    ttm  = mod_ttm,
    ttm_sep = mod_ttm_sep,
    ttm_cross = mod_ttm_cross
  )
  tab <- calc_loglik_tables(mods, cfg, S$X_te, config_canonical = config, perm = perm)
  # Remove bookkeeping column from main pipeline output
  if ("train_test_policy" %in% names(tab)) {
    tab$train_test_policy <- NULL
  }
  # Augment table with Kernel Smooth (marginal KDE) in nats
  LD_kern <- predict_kernel_marginal(mod_kern, S$X_te, type = "logdensity_by_dim")
  stopifnot(is.matrix(LD_kern), all(dim(LD_kern) == dim(S$X_te)))
  per_kern <- -colMeans(LD_kern)
  se_kern  <- apply(-LD_kern, 2, stderr)
  fmt <- function(m, se) sprintf("%.2f ± %.2f", round(m, 2), round(2 * se, 2))
  tab[["Kernel Smooth"]] <- c(
    fmt(per_kern, se_kern),
    fmt(sum(per_kern), stats::sd(rowSums(-LD_kern)) / sqrt(nrow(S$X_te)))
  )
  cat(sprintf("n=%d\n", n))
  print(tab)
  cat(sprintf("Permutation order %s (train/test only)\n", paste(perm, collapse = ",")))
  time_tab <- data.frame(
    model = c("True (marginal)", "True (Joint)", "Random Forest",
              "Kernel Smooth", "Marginal Map", "Separable Map", "Cross-term Map"),
    train_sec = c(t_true_tr, t_joint_tr, t_trtf_tr,
                  t_kern_tr, t_ttm_tr, t_sep_tr, t_ct_tr),
    test_sec = c(t_true_te, t_joint_te, t_trtf_te,
                 t_kern_te, t_ttm_te, t_sep_te, t_ct_te),
    stringsAsFactors = FALSE
  )
  time_tab$total_sec <- with(time_tab, train_sec + test_sec)
  stopifnot(all.equal(time_tab$total_sec,
                      time_tab$train_sec + time_tab$test_sec))
  print(time_tab)
  timing_table <<- time_tab
  results_table <<- tab
  stopifnot(identical(results_table, tab))
  return(tab)
}

if (sys.nframe() == 0L) {
  invisible(main())
}

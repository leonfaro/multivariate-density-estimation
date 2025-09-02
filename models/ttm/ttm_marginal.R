# Marginal Triangular Transport Map (TTM) — Minimal, monotone 1-D basis
# Implements fit_ttm, predict_ttm (with multiple types), and a tiny CLI for train/eval.

# try to source shared standardization util if available
if (!exists("standardize_train_only")) {
  p1 <- file.path("models", "ttm", "utils_standardize.R")
  if (file.exists(p1)) {
    source(p1)
  } else if (exists("root_path")) {
    p2 <- file.path(root_path, "models", "ttm", "utils_standardize.R")
    if (file.exists(p2)) source(p2)
  }
}

.std_stats <- function(X) {
  mu <- colMeans(X)
  sigma <- apply(X, 2, sd) + .Machine$double.eps
  list(mu = mu, sigma = sigma)
}

.standardize <- function(X, mu, sigma) {
  X1 <- sweep(X, 2, mu, "-")
  sweep(X1, 2, sigma, "/")
}

# Simple monotone f_k(x) = a_k + b_k x with b_k > 0
.fit_marginal_per_k <- function(x_std, eps = 1e-6) {
  # Optimize J(b,a) = 0.5 * sum (a + b x)^2 - n log b
  n <- length(x_std)
  xbar <- mean(x_std)
  v <- stats::var(x_std) + 1e-12
  b <- max(eps, 1 / sqrt(v))
  a <- -b * xbar
  c(a = a, b = b)
}

fit_ttm_marginal <- function(data, seed = 42) {
  set.seed(seed)
  S_in <- if (is.list(data) && !is.null(data$X_tr) && !is.null(data$X_te)) data else {
    stopifnot(is.matrix(data))
    if (!exists("split_data")) source("02_split.R")
    list(X_tr = split_data(data, seed)$X_tr, X_te = split_data(data, seed)$X_te)
  }
  X_tr <- S_in$X_tr; X_te <- S_in$X_te
  K <- ncol(X_tr)

  t_tr <- system.time({
    if (exists("standardize_train_only")) {
      st <- standardize_train_only(X_tr)
      mu <- st$mu; sigma <- st$sigma; X_tr_std <- st$Xs
    } else {
      st <- .std_stats(X_tr)
      mu <- st$mu; sigma <- st$sigma
      X_tr_std <- .standardize(X_tr, mu, sigma)
    }
    coeffs <- lapply(seq_len(K), function(k) .fit_marginal_per_k(X_tr_std[, k]))
    S <- list(algo = "marginal", mu = mu, sigma = sigma, coeffs = coeffs, order = seq_len(K))
    class(S) <- "ttm_marginal2"
  })[["elapsed"]]

  t_te <- system.time({ invisible(predict_ttm(S, X_te, type = "logdensity_by_dim")) })[["elapsed"]]

  list(S = S,
       NLL_train = mean(-predict_ttm(S, X_tr, type = "logdensity")),
       NLL_test  = mean(-predict_ttm(S, X_te,  type = "logdensity")),
       time_train = t_tr, time_pred = t_te)
}

# Dispatcher for different TTM algos
fit_ttm <- function(data, algo = c("marginal","separable","crossterm"), seed = 42, ...) {
  algo <- match.arg(algo)
  if (algo == "marginal") return(fit_ttm_marginal(data, seed = seed))
  if (algo == "separable") {
    src <- file.path("models", "ttm", "ttm_separable.R"); if (file.exists(src)) source(src)
    return(fit_ttm_separable(data, seed = seed, ...))
  }
  if (algo == "crossterm") {
    src <- file.path("models", "ttm", "ttm_crossterm.R"); if (file.exists(src)) source(src)
    return(fit_ttm_crossterm(data, seed = seed, ...))
  }
  stop(paste0("algo='", algo, "' not implemented in dispatcher"))
}

# Compatibility wrapper for existing scripts/tests
# Returns the same structure as fit_ttm_marginal (list with $S, times, NLLs)
trainMarginalMap <- function(S, seed = 42) {
  fit_ttm(S, algo = "marginal", seed = seed)
}

predict_ttm <- function(model, X, type = c("transform", "jac_diag", "logdensity_by_dim", "logdensity")) {
  type <- match.arg(type)
  # accept either a fitted list (with $S) or the model itself
  m <- if (is.list(model) && !is.null(model$S)) model$S else model
  if (!is.list(m) || is.null(m$mu) || is.null(m$sigma)) {
    stop("need TTM model or a fit with $S containing mu/sigma")
  }
  stopifnot(is.matrix(X))
  # Delegate to core
  if (!exists("ttm_forward")) {
    src <- if (exists("root_path")) file.path(root_path, "models", "ttm", "ttm_core.R") else file.path("models", "ttm", "ttm_core.R")
    if (file.exists(src)) source(src)
  }
  forward <- ttm_forward(m, X)
  if (type == "transform") return(forward$Z)
  if (type == "jac_diag") return(forward$J)
  LD <- ttm_ld_by_dim(m, X)
  stopifnot(is.matrix(LD), all(dim(LD) == dim(X)), all(is.finite(LD)))
  if (type == "logdensity_by_dim") return(LD)
  LDj <- rowSums(LD)
  stopifnot(all(is.finite(LDj)), max(abs(LDj - rowSums(LD))) <= 1e-10)
  LDj
}

# S3 predict wrapper
predict.ttm_marginal2 <- function(object, newdata, type = c("logdensity_by_dim", "logdensity"), ...) {
  predict_ttm(object, newdata, match.arg(type))
}

# Tiny CLI: Rscript R/ttm_marginal.R mode=train data=path.rds seed=42 out=mod.rds
train_val_test_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) return(invisible(NULL))
  kv <- strsplit(args, "=", fixed = TRUE)
  amap <- setNames(lapply(kv, function(x) paste(x[-1], collapse = "=")), vapply(kv, `[`, "", 1))
  mode <- amap[["mode"]] %||% "train_eval"
  seed <- as.integer(amap[["seed"]] %||% 42)
  data_path <- amap[["data"]]
  out <- amap[["out"]] %||% "ttm_marginal_model.rds"
  X <- NULL
  if (!is.null(data_path) && file.exists(data_path)) {
    obj <- readRDS(data_path)
    if (is.matrix(obj)) X <- obj else if (is.list(obj) && !is.null(obj$X)) X <- obj$X
  }
  if (is.null(X)) {
    N <- as.integer(amap[["N"]] %||% 400)
    K <- as.integer(amap[["K"]] %||% 3)
    set.seed(seed)
    X <- matrix(rnorm(N * K), ncol = K)
  }
  if (!exists("split_data")) source("02_split.R")
  S <- split_data(X, seed)
  if (mode %in% c("train", "train_eval", "fit")) {
    fit <- fit_ttm(S, seed = seed)
    saveRDS(fit$S, file = out)
    message(sprintf("[train] saved model to %s", out))
    if (mode == "train") return(invisible(TRUE))
  }
  if (mode %in% c("eval", "train_eval")) {
    M <- if (exists("fit") && is.list(fit)) fit$S else readRDS(out)
    nll <- mean(-predict_ttm(M, S$X_te, type = "logdensity"))
    cat(sprintf("NLL (nats): %.6f\n", nll))
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  # allow running as a script
  `%||%` <- function(a, b) if (is.null(a)) b else a
  train_val_test_cli()
}

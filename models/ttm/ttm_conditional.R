# Conditional sampling via truncated triangular inversion (composite reference)

# Helper: extract model (S) and optional training X
.extract_S_and_X <- function(model) {
  S <- if (is.list(model) && !is.null(model$S)) model$S else model
  X_tr <- NULL
  if (is.list(model) && !is.null(model$X_tr)) X_tr <- as.matrix(model$X_tr)
  if (is.list(S) && !is.null(S$X_tr)) X_tr <- as.matrix(S$X_tr)
  list(S = S, X_tr = X_tr)
}

# Core: invert k-th component for a single sample, given current x_prev (raw scale)
.invert_component_k <- function(S, k, z_target, x_prev_raw, bracket = c(-6, 6)) {
  mu <- S$mu; sigma <- S$sigma
  algo <- if (!is.null(S$algo)) S$algo else "marginal"
  K <- length(mu)
  # Fast path: marginal linear map
  if (identical(algo, "marginal")) {
    # Support both coeffs list (ttm_marginal2) and coeffA/coeffB legacy
    if (!is.null(S$coeffs)) {
      a <- as.numeric(S$coeffs[[k]]["a"])
      b <- as.numeric(S$coeffs[[k]]["b"])
    } else {
      b <- exp(S$coeffA[k]); a <- S$coeffB[k]
    }
    u_std <- (z_target - a) / b
    return(mu[k] + sigma[k] * u_std)
  }

  # Generic 1D monotone inversion via uniroot on standardized u
  f_eval <- function(u_std) {
    # Build a 1xK raw vector consistent with u_std at k and x_prev_raw for <k; others at training mean
    xr <- rep(0, K)
    # previous raw already supplied
    if (k > 1) xr[1:(k - 1)] <- as.numeric(x_prev_raw)
    xr[k] <- mu[k] + sigma[k] * u_std
    if (k < K) {
      # fill remaining with mu (u_std=0)
      xr[(k + 1):K] <- mu[(k + 1):K]
    }
    zk <- as.numeric(predict_ttm(S, matrix(xr, nrow = 1), type = "transform")[1, k])
    zk - z_target
  }
  lo <- bracket[1]; hi <- bracket[2]
  flo <- f_eval(lo); fhi <- f_eval(hi)
  # Expand bracket if no sign change (up to a few times)
  iters <- 0L
  while (flo * fhi > 0 && iters < 4L) {
    lo <- lo - 2; hi <- hi + 2
    flo <- f_eval(lo); fhi <- f_eval(hi)
    iters <- iters + 1L
  }
  # Final safeguard: if still no sign change, pick closer end
  if (flo * fhi > 0) {
    u_star <- if (abs(flo) < abs(fhi)) lo else hi
    return(mu[k] + sigma[k] * u_star)
  }
  rt <- uniroot(function(u) f_eval(u), lower = lo, upper = hi)
  mu[k] + sigma[k] * rt$root
}

#' Sample conditional using composite reference transport
#'
#' @param model TTM model or fit list with $S
#' @param x_cond numeric vector of conditioned values (raw scale)
#' @param idx_cond integer indices of conditioned coordinates (1-based)
#' @param n number of samples to draw
#' @param seed RNG seed for reproducibility
#' @return matrix of size n x |idx_free| in the original model order
sample_conditional_composite <- function(model, x_cond, idx_cond, n, seed = 42L) {
  stopifnot(is.numeric(x_cond), is.numeric(idx_cond), length(x_cond) == length(idx_cond))
  ex <- .extract_S_and_X(model); S <- ex$S; X_tr <- ex$X_tr
  if (!exists("predict_ttm")) {
    src <- file.path("models", "ttm", "ttm_marginal.R"); if (file.exists(src)) source(src)
  }
  if (is.null(S$mu) || is.null(S$sigma)) stop("Model must contain mu/sigma")
  K <- length(S$mu)
  idx_cond <- as.integer(idx_cond)
  if (length(idx_cond) == 0L) stop("idx_cond must be non-empty")
  if (any(idx_cond < 1 | idx_cond > K)) stop("idx_cond out of bounds")
  idx_free <- setdiff(seq_len(K), idx_cond)
  set.seed(as.integer(seed))
  # Choose source X rows
  if (is.null(X_tr)) stop("No training data found in model; attach as model$X_tr before calling")
  X_tr <- as.matrix(X_tr)
  stopifnot(ncol(X_tr) == K)
  sel <- sample.int(nrow(X_tr), size = n, replace = n > nrow(X_tr))
  X0 <- X_tr[sel, , drop = FALSE]
  # Composite reference Z_tilde
  Z_tilde <- predict_ttm(S, X0, type = "transform")
  out <- matrix(NA_real_, nrow = n, ncol = length(idx_free))
  colnames(out) <- paste0("x", idx_free)
  x_fix <- rep(NA_real_, K); x_fix[idx_cond] <- as.numeric(x_cond)
  for (i in seq_len(n)) {
    x_curr <- x_fix
    # Walk through triangular order k=1..K
    for (k in seq_len(K)) {
      if (!is.na(x_curr[k])) {
        # conditioned coordinate, keep fixed
        next
      }
      # build x_prev_raw (for j<k)
      if (k > 1) x_prev_raw <- x_curr[1:(k - 1)] else x_prev_raw <- numeric(0)
      z_target <- Z_tilde[i, k]
      xk_raw <- .invert_component_k(S, k, z_target, x_prev_raw)
      x_curr[k] <- xk_raw
    }
    out[i, ] <- x_curr[idx_free]
  }
  stopifnot(all(is.finite(out)))
  out
}


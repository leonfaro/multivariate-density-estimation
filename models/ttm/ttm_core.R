# Core helpers shared by TTM variants (marginal, separable, cross-term)

.ttm_std_apply <- function(mu, sigma, X) {
  X1 <- sweep(X, 2, mu, "-")
  sweep(X1, 2, sigma, "/")
}

# Minimal infix helper and config reader used across files
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a
.get_cross_maxit <- function() as.integer(getOption("cross.maxit", 200L))
.get_cross_quad_nodes <- function(fallback = NULL) {
  q <- getOption("cross.quad_nodes", NULL)
  if (!is.null(q)) return(as.integer(q))
  if (!is.null(fallback)) return(as.integer(fallback))
  as.integer(32L)
}
.get_cross_df_t <- function(fallback = NULL) {
  d <- getOption("cross.df_t", NULL)
  if (!is.null(d)) return(as.integer(d))
  if (!is.null(fallback)) return(as.integer(fallback))
  as.integer(8L)
}

.get_cross_lambda_non <- function(Q_used, N, fallback = NULL) {
  ln <- getOption("cross.lambda_non", NULL)
  if (!is.null(ln)) return(as.numeric(ln))
  if (!is.null(fallback) && is.finite(fallback)) return(as.numeric(fallback))
  0.05 * (as.numeric(Q_used) / as.numeric(N))
}

.get_cross_lambda_mon <- function(lam_non, Q_used, N, fallback = NULL) {
  lm <- getOption("cross.lambda_mon", NULL)
  if (!is.null(lm)) return(as.numeric(lm))
  if (!is.null(fallback) && is.finite(fallback)) return(as.numeric(fallback))
  0.5 * as.numeric(lam_non)
}

.use_cross_analytic_grad <- function() {
  isTRUE(getOption("cross.use_analytic_grad", TRUE))
}

.get_cross_workers <- function() {
  w <- getOption("cross.workers", NULL)
  if (!is.null(w)) return(as.integer(max(1L, w)))
  cw <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  as.integer(max(1L, cw - 1L))
}
.get_cross_deg_g <- function() as.integer(getOption("cross.deg_g", 3L))

# Determine algo tag from model object
.ttm_algo_of <- function(model) {
  if (!is.null(model$algo)) return(model$algo)
  cl <- class(model)
  if ("ttm_marginal2" %in% cl || "ttm_marginal" %in% cl) return("marginal")
  if ("ttm_separable" %in% cl) return("separable")
  if ("ttm_cross_term" %in% cl) return("crossterm")
  stop("Unknown TTM model class for core forward computation")
}

# ttm_forward: returns standardized forward map outputs and standardized Jacobian diagonal
# Z ∈ R^{N×K}, J ∈ R^{N×K} where J is ∂S/∂x_std (i.e., derivative w.r.t standardized input)
ttm_forward <- function(model, X) {
  stopifnot(is.matrix(X))
  mu <- model$mu; sigma <- model$sigma
  stopifnot(is.numeric(mu), is.numeric(sigma), length(mu) == ncol(X), length(sigma) == ncol(X))
  Xs <- .ttm_std_apply(mu, sigma, X)
  N <- nrow(Xs); K <- ncol(Xs)
  algo <- .ttm_algo_of(model)
  Z <- matrix(0, N, K)
  J <- matrix(0, N, K)
  if (identical(algo, "marginal")) {
    # Support both new (ttm_marginal2) and legacy (ttm_marginal) parameterizations
    if (!is.null(model$coeffs)) {
      for (k in seq_len(K)) {
        a <- as.numeric(model$coeffs[[k]]["a"])
        b <- as.numeric(model$coeffs[[k]]["b"])
        Z[, k] <- a + b * Xs[, k]
        J[, k] <- b
      }
    } else if (!is.null(model$coeffA) && !is.null(model$coeffB)) {
      b <- exp(model$coeffA)
      a <- model$coeffB
      for (k in seq_len(K)) {
        Z[, k] <- a[k] + b[k] * Xs[, k]
        J[, k] <- b[k]
      }
    } else {
      stop("Marginal model missing required coefficients")
    }
  } else if (identical(algo, "separable")) {
    # Build per-k using unified bases if available; fallback to legacy names
    if (!exists("build_g") || !exists("build_f") || !exists("d_build_f")) {
      # try load unified bases (new location)
      pth <- file.path("models", "ttm", "ttm_bases.R"); if (file.exists(pth)) source(pth)
    }
    use_new <- exists("build_g") && exists("build_f") && exists("d_build_f")
    deg <- if (!is.null(model$degree_g)) model$degree_g else 2L
    for (k in seq_len(K)) {
      x_prev <- if (k > 1) Xs[, 1:(k-1), drop = FALSE] else matrix(0, N, 0)
      xk <- Xs[, k]
      if (use_new) {
        if (k > 1) {
          P_non <- build_g(x_prev, deg)
        } else {
          # allow intercept-only for deg==0 to match marginal
          if (!is.null(model$degree_g) && model$degree_g == 0L && !is.null(model$coeffs[[k]]$c_non)) {
            P_non <- matrix(1, nrow = N, ncol = 1L)
          } else {
            P_non <- matrix(0, N, 0)
          }
        }
        P_mon <- build_f(xk)
        B <- d_build_f(xk)
        if (!is.null(model$coeffs[[k]]$c_mon) && length(model$coeffs[[k]]$c_mon) < ncol(P_mon)) {
          m_keep <- length(model$coeffs[[k]]$c_mon)
          P_mon <- P_mon[, seq_len(m_keep), drop = FALSE]
          B <- B[, seq_len(m_keep), drop = FALSE]
        }
      } else {
        if (!exists("basis_g") || !exists("basis_f") || !exists("dbasis_f")) {
          stop("Separable basis functions not found (build_* or basis_* set)")
        }
        P_non <- if (k > 1) basis_g(x_prev, deg) else matrix(0, N, 0)
        P_mon <- basis_f(xk)
        B <- dbasis_f(xk)
      }
      c_non <- model$coeffs[[k]]$c_non
      c_mon <- model$coeffs[[k]]$c_mon
      gk <- if (ncol(P_non) > 0) as.numeric(P_non %*% c_non) else rep(0, N)
      fk <- as.numeric(P_mon %*% c_mon)
      Z[, k] <- gk + fk
      J[, k] <- as.numeric(B %*% c_mon)
    }
  } else if (identical(algo, "crossterm")) {
    # Prefer unified bases build_g/build_h if present in model (new minimal crossterm)
    if (!is.null(model$spec_h) && !is.null(model$gl_nodes) && !is.null(model$gl_weights)) {
      nodes <- model$gl_nodes; weights <- model$gl_weights; Hmax <- if (!is.null(model$Hmax)) model$Hmax else 20
      for (k in seq_len(K)) {
        Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
        xk <- Xs[, k]
        Gx <- if (ncol(Xprev) > 0) build_g(Xprev, deg = model$deg_g) else matrix(1, N, 1L)
        alpha <- model$coeffs[[k]]$alpha
        beta  <- model$coeffs[[k]]$beta
        gk <- if (length(alpha) > 0) as.numeric(Gx %*% alpha) else rep(0, N)
        # Integral via GL on [0, xk]
        I <- rep(0, N)
        for (q in seq_along(nodes)) {
          tq <- xk * nodes[q]
          Hq <- build_h(tq, Xprev, model$spec_h)
          v <- as.numeric(Hq %*% beta)
          v <- pmax(pmin(v, Hmax), -Hmax)
          I <- I + weights[q] * exp(v)
        }
        I <- sign(xk) * abs(xk) * I
        Z[, k] <- gk + I
        # Jacobian diag: exp(h(xk, xprev))
        Hstar <- build_h(xk, Xprev, model$spec_h)
        hstar <- as.numeric(Hstar %*% beta)
        hstar <- pmax(pmin(hstar, Hmax), -Hmax)
        J[, k] <- exp(hstar)
      }
    } else {
      stop("Cross-term model not recognized by core; missing spec_h/gl_nodes")
    }
  } else stop("Unsupported algo")
  list(Z = Z, J = J)
}

# ttm_ld_by_dim: universal per-dim LD aggregator in log-space
ttm_ld_by_dim <- function(model, X) {
  out <- ttm_forward(model, X)
  Z <- out$Z; J <- out$J
  stopifnot(is.matrix(Z), is.matrix(J), all(dim(Z) == dim(J)))
  N <- nrow(Z); K <- ncol(Z)
  C <- -0.5 * log(2 * pi)
  LJ <- log(J) - matrix(log(model$sigma), nrow = N, ncol = K, byrow = TRUE)
  (-0.5) * (Z^2) + C + LJ
}


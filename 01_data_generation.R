#' Generate conditional samples
#'
#' This function implements conditional sampling via triangular transport mapping as described in `roadmap.md`.
#' It sequentially draws samples from the distributions specified in
#' `config`, potentially conditioning on previously generated columns.
#'

# 0.1  Bibliotheks-freie Pseudocode-Hilfsfunktionen -------------------

Generate_iid_from_config <- function(N , cfg){
  stopifnot(is.numeric(N), N > 0, is.list(cfg))
  K <- length(cfg)
  X <- matrix(NA_real_, nrow = N, ncol = K)
  for (i in seq_len(N)) {
    for (k in seq_len(K)) {
      c_k <- cfg[[k]]
      if (is.null(c_k$parm)) {
        args <- list()
      } else {
        if (k == 1) {
          prev <- data.frame()
        } else {
          prev <- as.data.frame(as.list(X[i, seq_len(k - 1)]))
          names(prev) <- paste0("X", seq_len(k - 1))
        }
        args <- c_k$parm(prev)
      }
      fun <- get(paste0("r", c_k$distr), mode = "function")
      if (c_k$distr == "gamma" &&
          all(c("shape1", "shape2") %in% names(args))) {
        args <- list(shape = args$shape1, scale = args$shape2)
      }
      args <- lapply(args, function(p) {
        if (!is.finite(p) || p <= 0) 1e-3 else p
      })
      X[i, k] <- do.call(fun, c(list(n = 1L), args))
    }
  }
  colnames(X) <- paste0("X", seq_len(K))
  X
}

Jacobian_Diagonal <- function(x, theta) {
  # Berechne ∂_{x_k} S_k  für k = 1..K  (theta enthaelt Map-Koeffizienten)
  diag(theta$L_inv)
}

LogDet_Jacobian <- function(x, theta) {
  # log|det ∇S| = Σ_k log( diagJ_k )
  sum(log(diag(theta$L_inv)))
}

S_forward <- function(x, theta) {
  # Untere-Dreieckige Abbildung:  S_k(x_1: k)
  as.numeric(theta$L_inv %*% (x - theta$mu))
}

R_inverse <- function(z, theta) {
  # Löse triangular:   x = R(z)
  as.numeric(theta$L %*% z + theta$mu)
}

Objective_J_N <- function(theta, X) {
  # J_N(theta) = -(1/N) Σ_n [ log η( S(x^(n);theta) ) + logDet_Jacobian(x^(n);theta ) ]
  stopifnot(is.list(theta), is.matrix(X))
  Z <- t(apply(X, 1L, S_forward, theta = theta))
  log_eta <- rowSums(dnorm(Z, log = TRUE))
  log_det <- rep(LogDet_Jacobian(X[1, ], theta), nrow(X))
  -(1 / nrow(X)) * sum(log_eta + log_det)
}

Conditional_Sample <- function(fix_idx, fix_val, theta_hat, m) {
  mu <- theta_hat$mu
  L <- theta_hat$L
  Sigma <- L %*% t(L)
  A <- setdiff(seq_along(mu), fix_idx)
  Sigma_BB <- Sigma[fix_idx, fix_idx, drop = FALSE]
  Sigma_AB <- Sigma[A, fix_idx, drop = FALSE]
  Sigma_AA <- Sigma[A, A, drop = FALSE]
  mu_B <- mu[fix_idx]
  mu_A <- mu[A]
  cond_mean <- c(mu_A + Sigma_AB %*% solve(Sigma_BB) %*% (fix_val - mu_B))
  cond_cov <- Sigma_AA - Sigma_AB %*% solve(Sigma_BB) %*% t(Sigma_AB)
  Z <- MASS::mvrnorm(n = m, mu = cond_mean, Sigma = cond_cov)
  if (m == 1) Z <- matrix(Z, nrow = 1)
  out <- matrix(NA_real_, nrow = m, ncol = length(mu))
  out[, fix_idx] <- matrix(rep(fix_val, each = m), nrow = m)
  out[, A] <- Z
  colnames(out) <- names(mu)
  out
}

# Kompatibilitäts-Funktion für bestehende Tests
gen_samples <- function(G) {
  Generate_iid_from_config(G$N, G$config)
}

#' Generate samples via Triangular Transport Map
#'
#' Implements the procedure described in `Theory.md` for a linear
#' triangular map.  First provisional data are drawn with the
#' sequential generator to fit the map parameters via BFGS on
#' `Objective_J_N`.  Final draws are obtained by applying the inverse
#' map to standard normal samples.
#'
#' @param config list of conditional specifications
#' @param N      integer sample size
#' @param seed   RNG seed ensuring reproducibility
#' @param fix_idx optional indices to condition on
#' @param fix_val numeric vector with fixed values
#' @param m       number of conditional samples if `fix_idx` is supplied
#' @return list with elements `X` and `theta_hat`; if `fix_idx` is given
#'   an additional matrix `X_cond` of conditional draws is returned
#' @export
TTM_generate <- function(config, N, seed, fix_idx = NULL, fix_val = NULL, m = 1L) {
  stopifnot(is.list(config), is.numeric(N), N > 0)
  set.seed(seed)
  K <- length(config)

  mu <- numeric(K)
  L_inv <- diag(1, K)


  N_fit <- ceiling(0.8 * N)
  X_fit <- Generate_iid_from_config(N_fit, config)

  idx_low  <- which(lower.tri(L_inv), arr.ind = TRUE)
  num_off  <- nrow(idx_low)
  idx_diag <- (K + num_off + 1):(K + num_off + K)
  par_init <- c(mu, L_inv[idx_low], rep(0, K))

  to_theta <- function(par) {
    mu_p      <- par[1:K]

    L_inv_p   <- diag(exp(par[idx_diag]), K)
    L_inv_p[idx_low] <- par[(K + 1):(K + num_off)]

    list(mu = mu_p, L_inv = L_inv_p)
  }

  fn <- function(par) Objective_J_N(to_theta(par), X_fit)

  gr <- function(par) {
    eps <- 1e-6
    sapply(seq_along(par), function(i) {
      p1 <- par; p1[i] <- p1[i] + eps
      p2 <- par; p2[i] <- p2[i] - eps
      (fn(p1) - fn(p2)) / (2 * eps)
    })
  }

  par <- par_init
  step_rel <- Inf
  grad_norm <- Inf
  iter <- 0L
  while (iter < 1000L && (grad_norm > 1e-6 || step_rel > 1e-8)) {
    iter <- iter + 1L
    opt <- optim(par, fn, gr, method = "BFGS", control = list(maxit = 1L, reltol = 0))
    par_new <- opt$par
    step_rel <- sqrt(sum((par_new - par)^2)) / (sqrt(sum(par^2)) + 1e-12)
    par <- par_new
    grad_curr <- gr(par)
    grad_norm <- sqrt(sum(grad_curr^2))
  }

  theta_hat <- to_theta(par)
  theta_hat$log_det_gradS <- sum(par[idx_diag])
  theta_hat$L  <- solve(theta_hat$L_inv)

  Z <- matrix(rnorm(K * N), nrow = N)
  X <- t(apply(Z, 1L, R_inverse, theta = theta_hat))

  colnames(X) <- paste0("X", seq_len(K))

  res <- list(X = X, theta_hat = theta_hat)
  if (!is.null(fix_idx)) {
    res$X_cond <- Conditional_Sample(fix_idx, fix_val, theta_hat = theta_hat, m)

  }
  res
}

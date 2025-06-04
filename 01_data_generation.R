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

Jacobian_Diagonal <- function(x , θ){
  # Berechne ∂_{x_k} S_k  für k = 1..K  (θ enthält Map-Koeffizienten)
  diag(θ$L_inv)
}

LogDet_Jacobian <- function(x , θ){
  # log|det ∇S| = Σ_k log( diagJ_k )
  sum(log(diag(θ$L_inv)))
}

S_forward <- function(x , θ){
  # Untere-Dreieckige Abbildung:  S_k(x_1: k)
  as.numeric(θ$L_inv %*% (x - θ$mu))
}

R_inverse <- function(z , θ){
  # Löse triangular:   x = R(z)
  as.numeric(θ$L %*% z + θ$mu)
}

Objective_J_N <- function(θ , X){
  # J_N(θ) = -(1/N) Σ_n [ log η( S(x^(n);θ) ) + logDet_Jacobian(x^(n);θ ) ]
  stopifnot(is.list(θ), is.matrix(X))
  Z <- t(apply(X, 1L, S_forward, θ = θ))
  log_eta <- rowSums(dnorm(Z, log = TRUE))
  log_det <- rep(LogDet_Jacobian(X[1, ], θ), nrow(X))
  -(1 / nrow(X)) * sum(log_eta + log_det)
}

Conditional_Sample <- function(fix_idx , fix_val , θ_hat , m){
  mu <- θ_hat$mu
  L <- θ_hat$L
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

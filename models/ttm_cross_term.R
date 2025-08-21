# Triangular Transport Map - Cross-Term Modul
# Basis-R Implementierung nach Vorgaben

# Hilfsfunktionen ------------------------------------------------------------

if (!exists(".standardizeData")) {
  .standardizeData <- function(X) {
    mu <- colMeans(X)
    sigma <- apply(X, 2, sd) + .Machine$double.eps
    X_tilde <- sweep(X, 2, mu, "-")
    X_tilde <- sweep(X_tilde, 2, sigma, "/")
    list(X = X_tilde, mu = mu, sigma = sigma)
  }
}

if (!exists(".standardize")) {
  .standardize <- function(S, X) {
    X <- sweep(X, 2, S$mu, "-")
    sweep(X, 2, S$sigma, "/")
  }
}

.basis_g_ct <- function(X, deg) {
  if (ncol(X) == 0L) {
    return(matrix(0, nrow = nrow(X), ncol = 0))
  }
  N <- nrow(X)
  out <- matrix(1, N, 1)
  for (j in seq_len(ncol(X))) {
    xj <- X[, j]
    for (d in seq_len(deg)) {
      out <- cbind(out, xj^d)
    }
  }
  out
}

.psi_basis_ct <- function(t, xprev, deg_t, deg_x, cross = TRUE,
                          deg_t_cross = 1, deg_x_cross = 1) {
  out <- numeric(0)
  if (deg_t > 0) {
    for (d in seq_len(deg_t)) {
      out <- c(out, t^d)
    }
  }
  if (cross && length(xprev) > 0) {
    for (j in seq_along(xprev)) {
      for (r in seq_len(deg_t_cross)) {
        for (s in seq_len(deg_x_cross)) {
          out <- c(out, t^r * xprev[j]^s)
        }
      }
    }
  }
  out
}

.dpsi_dt_ct <- function(t, xprev, deg_t, deg_x, cross = TRUE,
                        deg_t_cross = 1, deg_x_cross = 1) {
  out <- numeric(0)
  if (deg_t > 0) {
    for (d in seq_len(deg_t)) {
      out <- c(out, d * t^(max(d - 1, 0)))
    }
  }
  if (cross && length(xprev) > 0) {
    for (j in seq_along(xprev)) {
      for (r in seq_len(deg_t_cross)) {
        for (s in seq_len(deg_x_cross)) {
          out <- c(out, r * t^(max(r - 1, 0)) * xprev[j]^s)
        }
      }
    }
  }
  out
}

.build_Psi_q_ct <- function(xval, xp, nodes, nodes_pow, deg_t, deg_x,
                            deg_t_cross = 1, deg_x_cross = 1) {
  Q <- length(nodes)
  m_beta <- deg_t + length(xp) * deg_t_cross * deg_x_cross
  Psi_q <- matrix(0, Q, m_beta)
  if (deg_t > 0) {
    x_pow <- xval^(seq_len(deg_t))
    Psi_q[, seq_len(deg_t)] <- sweep(nodes_pow[, seq_len(deg_t), drop = FALSE], 2, x_pow, "*")
  }
  if (length(xp) > 0) {
    col <- deg_t
    x_pow_cross <- xval^(seq_len(deg_t_cross))
    xp_pows <- lapply(xp, function(xj) xj^(seq_len(deg_x_cross)))
    for (j in seq_along(xp)) {
      for (r in seq_len(deg_t_cross)) {
        for (s in seq_len(deg_x_cross)) {
          col <- col + 1
          Psi_q[, col] <- nodes_pow[, r] * x_pow_cross[r] * xp_pows[[j]][s]
        }
      }
    }
  }
  Psi_q
}
.gauss_legendre_01_ct <- function(n) {
  if (n <= 0 || n != as.integer(n)) {
    stop("n must be positive integer")
  }
  if (n == 1) {
    return(list(nodes = 0.5, weights = 1))
  }
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

.safe_mclapply_ct <- function(X, FUN, mc.cores, ...) {
  parallel::mclapply(X, function(ix) {
    tryCatch(FUN(ix), error = function(e) e)
  }, mc.cores = mc.cores, mc.set.seed = FALSE, mc.preschedule = TRUE, ...)
}

.is_coeffs_ok <- function(x) is.list(x) && is.numeric(x$alpha) && is.numeric(x$beta)
.is_predchunk_ok <- function(x) is.list(x) && is.numeric(x$Z_col) && is.numeric(x$LJ_col)

.zero_coeffs_ct <- function(k, X_tr_std, degree_g, degree_t,
                           degree_t_cross, degree_x_cross) {
  N <- nrow(X_tr_std)
  Xprev <- if (k > 1) X_tr_std[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
  m_alpha <- ncol(.basis_g_ct(Xprev, degree_g))
  xprev_first <- if (k > 1) Xprev[1, , drop = TRUE] else numeric(0)
        m_beta <- length(.psi_basis_ct(0, xprev_first, degree_t, degree_g, TRUE, degree_t_cross, degree_x_cross))
  list(alpha = if (m_alpha > 0) rep(0, m_alpha) else numeric(0),
       beta  = rep(0, m_beta),
       convergence = NA_real_)
}

.zero_predchunk_ct <- function(N, sigma_k) {
  list(Z_col = rep(0, N), LJ_col = rep(-log(sigma_k), N))
}

.NLL_set_ct <- function(S, X) {
  mean(-rowSums(predict(S, X, "logdensity_by_dim")))
}

.SE_set_ct <- function(S, X) {
  v <- rowSums(-predict(S, X, "logdensity_by_dim"))
  stats::sd(v) / sqrt(length(v))
}

.logJacDiag_ct <- function(S, x) {
  Xs <- .standardize(S, matrix(x, nrow = 1))
  K <- length(x)
  out <- numeric(K)
  for (k in seq_len(K)) {
    xprev <- if (k > 1) Xs[1, 1:(k - 1)] else numeric(0)
    psi <- .psi_basis_ct(Xs[1, k], xprev, S$degree_t, S$degree_g, TRUE, S$degree_t_cross, S$degree_x_cross)
    beta <- S$coeffs[[k]]$beta
    out[k] <- sum(beta * psi) - log(S$sigma[k])
  }
  out
}

.forwardKLLoss_ct <- function(S, X) {
  X <- as.matrix(X)
  K <- ncol(X)
  LD <- predict(S, X, "logdensity_by_dim")
  mean(-rowSums(LD) - 0.5 * K * log(2 * pi))
}

# Exportierte Funktionen -----------------------------------------------------

trainCrossTermMap <- function(X_or_path, degree_g = 2, degree_t = 2, degree_t_cross = 1, degree_x_cross = 1,
                              lambda = 1e-3, batch_n = NULL, Q = NULL,
                              eps = 1e-6, clip = Inf,
                              alpha_init_list = NULL, warmstart_from_separable = FALSE,
                              sep_degree_g = NULL, sep_lambda = 1e-3) {
  set.seed(42)
  S_in <- if (is.character(X_or_path)) readRDS(X_or_path) else X_or_path
  stopifnot(is.list(S_in))
  X_tr <- S_in$X_tr
  X_val <- S_in$X_val
  X_te  <- S_in$X_te

  if (is.null(alpha_init_list) && warmstart_from_separable) {
    if (!exists("trainSeparableMap")) {
      stop("trainSeparableMap not found for warm start")
    }
    sep_deg <- if (is.null(sep_degree_g)) degree_g else sep_degree_g
    fit_sep <- trainSeparableMap(S_in, degree_g = sep_deg, lambda = sep_lambda)
    alpha_init_list <- lapply(fit_sep$S$coeffs, `[[`, "c_non")
  }

  time_train <- system.time({
    std <- .standardizeData(X_tr)
    X_tr_std <- std$X
    mu <- std$mu
    sigma <- std$sigma
    N <- nrow(X_tr_std)
    K <- ncol(X_tr_std)
    if (!is.null(alpha_init_list)) {
      stopifnot(length(alpha_init_list) == K)
    }
      degree_t_max <- max(degree_t, degree_t_cross)
      Q_use <- if (is.null(Q)) min(12, 4 + 2 * degree_t_max) else Q
      batch_use <- if (is.null(batch_n)) min(N, max(256L, floor(65536 / max(1L, Q_use)))) else min(N, batch_n)
    quad <- .gauss_legendre_01_ct(Q_use)
    nodes <- quad$nodes
    weights <- quad$weights
      nodes_pow <- outer(nodes, seq_len(degree_t_max), `^`)

    fit_k <- function(k) {
      Xprev <- if (k > 1) X_tr_std[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
      xk <- X_tr_std[, k]
      Phi <- .basis_g_ct(Xprev, degree_g)
      m_alpha <- ncol(Phi)
      xprev_first <- if (k > 1) Xprev[1, , drop = TRUE] else numeric(0)
        m_beta <- length(.psi_basis_ct(0, xprev_first, degree_t, degree_g, TRUE, degree_t_cross, degree_x_cross))
      alpha_start <- if (m_alpha > 0) {
        if (is.null(alpha_init_list)) {
          rep(0, m_alpha)
        } else {
          ai <- alpha_init_list[[k]]
          if (length(ai) != m_alpha) rep(0, m_alpha) else ai
        }
      } else {
        numeric(0)
      }

      loss_grad <- function(alpha, beta) {
        S_sq_sum <- 0
        term_sum <- 0
        grad_alpha <- if (m_alpha > 0) rep(0, m_alpha) else numeric(0)
        grad_beta <- rep(0, m_beta)
        for (i0 in seq(1, N, by = batch_use)) {
          idx <- i0:min(i0 + batch_use - 1, N)
          Phi_blk <- if (m_alpha > 0) Phi[idx, , drop = FALSE] else NULL
          Xprev_blk <- if (k > 1) Xprev[idx, , drop = FALSE] else matrix(0, length(idx), 0)
          xk_blk <- xk[idx]
          for (b in seq_along(idx)) {
            xp <- if (k > 1) Xprev_blk[b, ] else numeric(0)
            xval <- xk_blk[b]
            Psi_q <- .build_Psi_q_ct(xval, xp, nodes, nodes_pow, degree_t, degree_g, degree_t_cross, degree_x_cross)
            V <- as.vector(Psi_q %*% beta)
            b_vec <- log(weights) + V
            b_max <- max(b_vec)
            r <- exp(b_vec - b_max)
            s <- exp(b_max) * sum(r)
            I_i <- xval * s
            soft <- r / sum(r)
            dI_i <- I_i * as.vector(t(Psi_q) %*% soft)
            psi_x <- .psi_basis_ct(xval, xp, degree_t, degree_g, TRUE, degree_t_cross, degree_x_cross)
            S_i <- if (m_alpha > 0) sum(Phi_blk[b, ] * alpha) + I_i else I_i
            S_sq_sum <- S_sq_sum + S_i^2
            if (m_alpha > 0) grad_alpha <- grad_alpha + S_i * Phi_blk[b, ]
            grad_beta <- grad_beta + S_i * dI_i - psi_x
            term_sum <- term_sum + sum(psi_x * beta)
          }
        }
        loss <- 0.5 * S_sq_sum / N - term_sum / N +
          0.5 * lambda * (sum(alpha^2) + sum(beta^2))
        grad_alpha_out <- if (m_alpha > 0) grad_alpha / N + lambda * alpha else numeric(0)
        grad_beta_out <- grad_beta / N + lambda * beta
        list(loss = loss, grad = c(grad_alpha_out, grad_beta_out))
      }

      cache <- new.env(parent = emptyenv())
      fn <- function(theta) {
        alpha <- if (m_alpha > 0) theta[seq_len(m_alpha)] else numeric(0)
        beta <- theta[(m_alpha + 1):length(theta)]
        res <- loss_grad(alpha, beta)
        cache$grad <- res$grad
        res$loss
      }
      gr <- function(theta) {
        if (!is.null(cache$grad)) {
          g <- cache$grad
          cache$grad <- NULL
          return(g)
        }
        alpha <- if (m_alpha > 0) theta[seq_len(m_alpha)] else numeric(0)
        beta <- theta[(m_alpha + 1):length(theta)]
        loss_grad(alpha, beta)$grad
      }

      theta0 <- c(alpha_start, rep(0, m_beta))
      opt <- optim(theta0, fn, gr, method = "L-BFGS-B", lower = -Inf, upper = Inf)
      alpha_hat <- if (m_alpha > 0) opt$par[seq_len(m_alpha)] else numeric(0)
      beta_hat <- opt$par[(m_alpha + 1):length(opt$par)]
      list(alpha = alpha_hat, beta = beta_hat, convergence = opt$convergence)
    }

    res <- .safe_mclapply_ct(seq_len(K), fit_k,
      mc.cores = min(getOption("mc.cores", 10L), K))
    for (k in seq_len(K)) {
      if (!.is_coeffs_ok(res[[k]])) {
        r2 <- tryCatch(fit_k(k), error = function(e) e)
        if (!.is_coeffs_ok(r2)) {
            r2 <- .zero_coeffs_ct(k, X_tr_std, degree_g, degree_t, degree_t_cross, degree_x_cross)
        }
        res[[k]] <- r2
      }
    }
    coeffs <- res
    stopifnot(length(coeffs) == K,
              all(vapply(coeffs, .is_coeffs_ok, logical(1))) )

    S_map <- list(
      mu = mu,
      sigma = sigma,
      degree_g = degree_g,
      degree_t = degree_t,
      degree_t_cross = degree_t_cross,
      degree_x_cross = degree_x_cross,
      coeffs = coeffs,
      Q = Q_use,
      nodes = nodes,
      weights = weights,
      nodes_pow = nodes_pow,
      quad_nodes_ct = nodes,
      quad_weights_ct = weights,
      quad_nodes_pow_ct = nodes_pow,
      clip = clip,
      order = seq_len(K)
    )
    class(S_map) <- "ttm_cross_term"
  })[["elapsed"]]

  time_pred <- system.time({
    predict(S_map, X_te, "logdensity_by_dim")
  })[["elapsed"]]

list(
    S = S_map,
    NLL_train = .NLL_set_ct(S_map, X_tr),
    NLL_val = .NLL_set_ct(S_map, X_val),
    NLL_test = .NLL_set_ct(S_map, X_te),
    stderr_test = .SE_set_ct(S_map, X_te),
    time_train = time_train,
    time_pred = time_pred
  )
}

predict.ttm_cross_term <- function(object, newdata,
                                   type = c("logdensity_by_dim", "logdensity"),
                                   batch_n = NULL) {
  type <- tryCatch(match.arg(type), error = function(e) stop("unknown type"))
  Xs <- .standardize(object, newdata)
  N <- nrow(Xs)
  K <- ncol(Xs)
  batch_use <- if (is.null(batch_n)) min(N, max(256L, floor(65536 / max(1L, object$Q)))) else min(N, batch_n)
  nodes <- object$quad_nodes_ct
  weights <- object$quad_weights_ct
  nodes_pow <- object$quad_nodes_pow_ct
  C <- -0.5 * log(2 * pi)

  chunk_fun <- function(k) {
    Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
    xk <- Xs[, k]
    Phi <- .basis_g_ct(Xprev, object$degree_g)
    m_alpha <- ncol(Phi)
    alpha <- object$coeffs[[k]]$alpha
    beta <- object$coeffs[[k]]$beta
    xprev_first <- if (k > 1) Xprev[1, , drop = TRUE] else numeric(0)
      psi_len <- length(.psi_basis_ct(Xs[1, k], xprev_first,
                                      object$degree_t, object$degree_g, TRUE, object$degree_t_cross, object$degree_x_cross))
    stopifnot(length(beta) == psi_len)
    Z_col <- numeric(N)
    LJ_col <- numeric(N)
    for (i0 in seq(1, N, by = batch_use)) {
      idx <- i0:min(i0 + batch_use - 1, N)
      Phi_blk <- if (m_alpha > 0) Phi[idx, , drop = FALSE] else NULL
      Xprev_blk <- if (k > 1) Xprev[idx, , drop = FALSE] else matrix(0, length(idx), 0)
      xk_blk <- xk[idx]
      for (b in seq_along(idx)) {
        xp <- if (k > 1) Xprev_blk[b, ] else numeric(0)
        xval <- xk_blk[b]
        Psi_q <- .build_Psi_q_ct(xval, xp, nodes, nodes_pow, object$degree_t, object$degree_g, object$degree_t_cross, object$degree_x_cross)
        V <- as.vector(Psi_q %*% beta)
        b_vec <- log(weights) + V
        b_max <- max(b_vec)
        r <- exp(b_vec - b_max)
        s <- exp(b_max) * sum(r)
        I_i <- xval * s
        psi_x <- .psi_basis_ct(xval, xp, object$degree_t, object$degree_g, TRUE, object$degree_t_cross, object$degree_x_cross)
        Z_col[idx[b]] <- if (m_alpha > 0) sum(Phi_blk[b, ] * alpha) + I_i else I_i
        LJ_col[idx[b]] <- sum(beta * psi_x) - log(object$sigma[k])
      }
    }
    list(Z_col = Z_col, LJ_col = LJ_col)
  }

  res <- .safe_mclapply_ct(seq_len(K), chunk_fun,
    mc.cores = min(getOption("mc.cores", 10L), K))
  for (k in seq_len(K)) {
    if (!.is_predchunk_ok(res[[k]])) {
      r2 <- tryCatch(chunk_fun(k), error = function(e) e)
      if (!.is_predchunk_ok(r2)) {
        r2 <- .zero_predchunk_ct(N, object$sigma[k])
        r2$Z_col <- Xs[, k]
      }
      res[[k]] <- r2
    }
  }

  Z <- do.call(cbind, lapply(res, `[[`, "Z_col"))
  LJ <- do.call(cbind, lapply(res, `[[`, "LJ_col"))
  LD <- (-0.5) * (Z^2) + C + LJ
  if (type == "logdensity_by_dim") {
    LD
  } else {
    rowSums(LD)
  }
}


# Ende -----------------------------------------------------------------------

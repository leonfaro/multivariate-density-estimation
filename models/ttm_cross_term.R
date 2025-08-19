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

.psi_basis_ct <- function(t, xprev, deg_t, deg_x, cross = TRUE) {
  out <- numeric(0)
  for (d in seq_len(deg_t)) {
    out <- c(out, t^d)
  }
  if (length(xprev) > 0) {
    for (j in seq_along(xprev)) {
      for (d in seq_len(deg_x)) {
        out <- c(out, xprev[j]^d)
      }
    }
    if (cross) {
      for (j in seq_along(xprev)) {
        out <- c(out, t * xprev[j])
      }
    }
  }
  out
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
    psi <- .psi_basis_ct(Xs[1, k], xprev, S$degree_t, S$degree_g, TRUE)
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

trainCrossTermMap <- function(X_or_path, degree_g = 2, degree_t = 2,
                              lambda = 1e-3, Q = 16, eps = 1e-6,
                              clip = 20) {
  set.seed(42)
  S_in <- if (is.character(X_or_path)) readRDS(X_or_path) else X_or_path
  stopifnot(is.list(S_in))
  X_tr <- S_in$X_tr
  X_val <- S_in$X_val
  X_te  <- S_in$X_te

  time_train <- system.time({
    std <- .standardizeData(X_tr)
    X_tr_std <- std$X
    mu <- std$mu
    sigma <- std$sigma
    N <- nrow(X_tr_std)
    K <- ncol(X_tr_std)
    coeffs <- vector("list", K)

    quad <- .gauss_legendre_01_ct(Q)
    nodes <- quad$nodes
    weights <- quad$weights

    for (k in seq_len(K)) {
      Xprev <- if (k > 1) X_tr_std[, 1:(k - 1), drop = FALSE] else
        matrix(0, N, 0)
      xk <- X_tr_std[, k]
      Phi <- .basis_g_ct(Xprev, degree_g)
      m_alpha <- ncol(Phi)
      xprev_first <- if (k > 1) Xprev[1, , drop = TRUE] else numeric(0)
      m_beta <- length(.psi_basis_ct(0, xprev_first, degree_t, degree_g, TRUE))
      Psi_x <- matrix(0, N, m_beta)
      for (i in seq_len(N)) {
        xp <- if (k > 1) Xprev[i, , drop = TRUE] else numeric(0)
        Psi_x[i, ] <- .psi_basis_ct(xk[i], xp, degree_t, degree_g, TRUE)
      }

      calc_I <- function(beta) {
        I <- numeric(N)
        dI <- matrix(0, N, m_beta)
        for (i in seq_len(N)) {
          xp <- if (k > 1) Xprev[i, , drop = TRUE] else numeric(0)
          xval <- xk[i]
          t_vec <- nodes * xval
          Psi_q <- t(vapply(t_vec, function(tt)
            .psi_basis_ct(tt, xp, degree_t, degree_g, TRUE),
            numeric(m_beta)))
          v <- Psi_q %*% beta
          v <- pmin(pmax(v, -clip), clip)
          e <- exp(v)
          I[i] <- xval * sum(weights * e)
          dI[i, ] <- xval * as.vector(t(Psi_q) %*% (weights * e))
        }
        list(I = I, dI = dI)
      }

      fn <- function(theta) {
        alpha <- if (m_alpha > 0) theta[seq_len(m_alpha)] else numeric(0)
        beta <- theta[(m_alpha + 1):length(theta)]
        calc <- calc_I(beta)
        S_vec <- calc$I
        if (m_alpha > 0) {
          S_vec <- S_vec + as.vector(Phi %*% alpha)
        }
        term <- Psi_x %*% beta
        mean(0.5 * S_vec^2 - term) +
          0.5 * lambda * (sum(alpha^2) + sum(beta^2))
      }

      gr <- function(theta) {
        alpha <- if (m_alpha > 0) theta[seq_len(m_alpha)] else numeric(0)
        beta <- theta[(m_alpha + 1):length(theta)]
        calc <- calc_I(beta)
        S_vec <- calc$I
        if (m_alpha > 0) {
          S_vec <- S_vec + as.vector(Phi %*% alpha)
        }
        grad_alpha <- if (m_alpha > 0) colMeans(S_vec * Phi) + lambda * alpha else numeric(0)
        grad_beta <- colMeans(S_vec * calc$dI) -
          colMeans(Psi_x) + lambda * beta
        c(grad_alpha, grad_beta)
      }

      theta0 <- rep(0, m_alpha + m_beta)
      opt <- optim(theta0, fn, gr, method = "L-BFGS-B",
                   control = list(fnscale = 1), lower = -Inf, upper = Inf)
      alpha_hat <- if (m_alpha > 0) opt$par[seq_len(m_alpha)] else numeric(0)
      beta_hat <- opt$par[(m_alpha + 1):length(opt$par)]
      coeffs[[k]] <- list(alpha = alpha_hat, beta = beta_hat)
    }

    S_map <- list(
      mu = mu,
      sigma = sigma,
      degree_g = degree_g,
      degree_t = degree_t,
      coeffs = coeffs,
      Q = Q,
      quad_nodes_ct = nodes,
      quad_weights_ct = weights,
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
                                   type = c("logdensity_by_dim", "logdensity")) {
  type <- tryCatch(match.arg(type), error = function(e) stop("unknown type"))
  Xs <- .standardize(object, newdata)
  N <- nrow(Xs)
  K <- ncol(Xs)
  Z <- matrix(0, N, K)
  LJ <- matrix(0, N, K)
  nodes <- object$quad_nodes_ct
  weights <- object$quad_weights_ct
  C <- -0.5 * log(2 * pi)

  for (k in seq_len(K)) {
    Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
    xk <- Xs[, k]
    Phi <- .basis_g_ct(Xprev, object$degree_g)
    m_alpha <- ncol(Phi)
    alpha <- object$coeffs[[k]]$alpha
    beta <- object$coeffs[[k]]$beta
    xprev_first <- if (k > 1) Xs[1, 1:(k - 1), drop = TRUE] else numeric(0)
    psi_len <- length(.psi_basis_ct(Xs[1, k], xprev_first,
                                    object$degree_t, object$degree_g, TRUE))
    stopifnot(length(beta) == psi_len)
    m_beta <- length(beta)
    I <- numeric(N)
    psi_x <- matrix(0, N, m_beta)
    for (i in seq_len(N)) {
      xp <- if (k > 1) Xprev[i, , drop = TRUE] else numeric(0)
      xval <- xk[i]
      t_vec <- nodes * xval
      Psi_q <- t(vapply(t_vec, function(tt)
        .psi_basis_ct(tt, xp, object$degree_t, object$degree_g, TRUE),
        numeric(m_beta)))
      v <- Psi_q %*% beta
      v <- pmin(pmax(v, -object$clip), object$clip)
      e <- exp(v)
      I[i] <- xval * sum(weights * e)
      psi_x[i, ] <- .psi_basis_ct(xval, xp, object$degree_t, object$degree_g, TRUE)
    }
    if (m_alpha > 0) {
      Z[, k] <- as.vector(Phi %*% alpha) + I
    } else {
      Z[, k] <- I
    }
    LJ[, k] <- psi_x %*% beta - log(object$sigma[k])
  }

  LD <- (-0.5) * (Z^2) + C + LJ
  if (type == "logdensity_by_dim") {
    LD
  } else {
    rowSums(LD)
  }
}

# Ende -----------------------------------------------------------------------

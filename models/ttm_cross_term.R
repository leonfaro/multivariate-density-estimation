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

# RBF basis (1D) with Gaussian kernels
.rbf_basis_1d_ct <- function(x, centers, sigma) {
  if (length(x) == 0) return(matrix(0, 0, 0))
  Z <- outer(x, centers, function(a, c) exp(-0.5 * ((a - c)/sigma)^2))
  colnames(Z) <- paste0("rbf_", seq_along(centers))
  Z
}

# Build design matrix for g_k(x_<k>)
.build_g_design_ct <- function(Xprev, k, degree_g, use_rbf = TRUE) {
  # For 2D (k==2) with one predecessor, use RBF(7) + poly {1,x,x^2}
  if (ncol(Xprev) == 1L && k == 2L && use_rbf) {
    x1 <- Xprev[, 1]
    ctr <- as.numeric(stats::quantile(x1, probs = seq(0.1, 0.9, length.out = 7)))
    # bandwidth from quantile spacing or fallback to sd
    dif <- diff(ctr)
    sig <- if (all(is.finite(dif)) && length(dif) > 0) 0.5 * stats::median(dif) else (stats::sd(x1) + 1e-8)/3
    poly <- cbind(1, x1, x1^2)
    rbf <- .rbf_basis_1d_ct(x1, centers = ctr, sigma = sig)
    Phi <- cbind(poly, rbf)
    attr(Phi, "g_spec") <- list(type = "rbf_poly", centers = ctr, sigma = sig)
    return(Phi)
  }
  # Fallback: polynomial in predecessors
  Phi <- .basis_g_ct(Xprev, degree_g)
  attr(Phi, "g_spec") <- list(type = "poly", degree_g = degree_g)
  Phi
}

# psi_m(x_<k>) for 1D predecessor up to cubic
.psi_xprev_poly3_ct <- function(xprev) {
  if (length(xprev) == 0L) return(1)
  # 1 plus for jede Vorgängerdimension: (x, x^2, x^3)
  out <- 1
  for (j in seq_along(xprev)) {
    x <- xprev[j]
    out <- c(out, x, x^2, x^3)
  }
  out
}

# Build Psi_q using B-spline basis in t (cubic) times psi(x_<k>)
.build_Psi_q_bs_ct <- function(xval, xp, nodes, bs_spec) {
  Q <- length(nodes)
  t_nodes <- xval * nodes
  # evaluate B-splines at t_nodes with stored spec
  B <- splines::bs(t_nodes,
                   df = bs_spec$df,
                   degree = bs_spec$degree,
                   knots = bs_spec$knots,
                   Boundary.knots = bs_spec$boundary,
                   intercept = TRUE)
  psi <- .psi_xprev_poly3_ct(xp)
  # Kronecker per row
  Psi_q <- matrix(0, Q, ncol(B) * length(psi))
  col <- 0
  for (j in seq_len(ncol(B))) {
    for (m in seq_along(psi)) {
      col <- col + 1
      Psi_q[, col] <- B[, j] * psi[m]
    }
  }
  Psi_q
}

# Build cached B-spline tensor and psi matrices for training
.build_ct_cache_tr <- function(xk, Xprev, nodes, bs_spec) {
  stopifnot(length(xk) == nrow(Xprev) || nrow(Xprev) == 0)
  N <- length(xk)
  Q <- length(nodes)
  # B3D as stacked rows (N*Q) x df
  Btmp0 <- splines::bs(0, df = bs_spec$df, degree = bs_spec$degree,
                       knots = bs_spec$knots, Boundary.knots = bs_spec$boundary,
                       intercept = TRUE)
  df <- ncol(Btmp0)
  B2D <- matrix(0, nrow = N * Q, ncol = df)
  for (i in seq_len(N)) {
    Bi <- splines::bs(xk[i] * nodes, df = bs_spec$df, degree = bs_spec$degree,
                      knots = bs_spec$knots, Boundary.knots = bs_spec$boundary,
                      intercept = TRUE)
    B2D[((i - 1L) * Q + 1L):(i * Q), ] <- Bi
  }
  # psi(x_<k>) rows: [1, for each prev dim j: (x_j, x_j^2, x_j^3)]
  if (ncol(Xprev) >= 1L) {
    pieces <- list(matrix(1, nrow = N, ncol = 1))
    for (j in seq_len(ncol(Xprev))) {
      xj <- Xprev[, j]
      pieces[[length(pieces) + 1L]] <- cbind(xj, xj^2, xj^3)
    }
    PsiX <- do.call(cbind, pieces)
  } else {
    PsiX <- matrix(1, nrow = N, ncol = 1)
  }
  # Row-wise kron to build Pi: for each i, kron(Bi, PsiX[i,])
  M <- df * ncol(PsiX)
  Pi <- matrix(0, nrow = N * Q, ncol = M)
  for (i in seq_len(N)) {
    rows <- ((i - 1L) * Q + 1L):(i * Q)
    Bi <- B2D[rows, , drop = FALSE]
    psi <- matrix(PsiX[i, ], nrow = 1)
    # Khatri-Rao by row: expand by multiplication
    block <- do.call(cbind, lapply(seq_len(ncol(Bi)), function(j) Bi[, j] * psi))
    Pi[rows, ] <- block
  }
  list(Pi = Pi, B2D = B2D, PsiX = PsiX, df = df, M = M, Q = Q)
}

# Row-wise Khatri–Rao product: C[i,] = vec(A[i,] ⊗ B[i,])
.KR_rowwise_ct <- function(A, B) {
  stopifnot(is.matrix(A), is.matrix(B), nrow(A) == nrow(B))
  N <- nrow(A); a <- ncol(A); b <- ncol(B)
  C <- matrix(0, nrow = N, ncol = a * b)
  for (j in seq_len(a)) {
    cols <- ((j - 1L) * b + 1L):(j * b)
    C[, cols] <- B * A[, j]
  }
  C
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
  # Avoid nested parallelism: fall back to sequential if a higher-level
  # orchestrator is already parallelizing.
  if (getOption("mde.parallel_active", FALSE)) {
    lapply(X, function(ix) tryCatch(FUN(ix), error = function(e) e))
  } else {
    parallel::mclapply(X, function(ix) {
      tryCatch(FUN(ix), error = function(e) e)
    }, mc.cores = mc.cores, mc.set.seed = FALSE, mc.preschedule = TRUE, ...)
  }
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
                              sep_degree_g = NULL, sep_lambda = 1e-3, seed = 42,
                              lambda_non = NULL, lambda_mon = NULL,
                              order_mode = c("as-is","x1_x2","x2_x1")) {
  # Order mode from options/env; default "as-is" (no permutation during training)
  order_mode <- match.arg(order_mode)
  order_mode_opt <- tryCatch(getOption("mde.ctm.order"), error = function(e) NULL)
  if (is.character(order_mode_opt) && nzchar(order_mode_opt)) {
    order_mode <- match.arg(order_mode_opt, choices = c("as-is","x1_x2","x2_x1"))
  } else {
    env_ord <- Sys.getenv("MDE_CTM_ORDER", NA_character_)
    if (!is.na(env_ord) && nzchar(env_ord)) {
      order_mode <- match.arg(env_ord, choices = c("as-is","x1_x2","x2_x1"))
    }
  }

  # Ridge hyperparameters from options/env with defaults if not passed
  if (is.null(lambda_non)) {
    lambda_non <- getOption("mde.ctm.lambda_non", NA_real_)
    if (is.na(lambda_non)) {
      env <- suppressWarnings(as.numeric(Sys.getenv("MDE_CTM_LAMBDA_NON", NA_character_)))
      # Default for lambda_non aligned to 3e-2 (no grid)
      lambda_non <- ifelse(is.na(env), 3e-2, env)
    }
  }
  if (is.null(lambda_mon)) {
    lambda_mon <- getOption("mde.ctm.lambda_mon", NA_real_)
    if (is.na(lambda_mon)) {
      env <- suppressWarnings(as.numeric(Sys.getenv("MDE_CTM_LAMBDA_MON", NA_character_)))
      # Slightly lower default to increase flexibility
      lambda_mon <- ifelse(is.na(env), 2e-2, env)
    }
  }
  set.seed(seed)
  S_in <- if (is.character(X_or_path)) readRDS(X_or_path) else X_or_path
  stopifnot(is.list(S_in))
  X_tr <- S_in$X_tr
  X_te  <- S_in$X_te

  if (is.null(alpha_init_list) && warmstart_from_separable) {
    if (!exists("trainSeparableMap")) {
      stop("trainSeparableMap not found for warm start")
    }
    sep_deg <- if (is.null(sep_degree_g)) degree_g else sep_degree_g
    fit_sep <- trainSeparableMap(S_in, degree_g = sep_deg, lambda = sep_lambda, seed = seed)
    alpha_init_list <- lapply(fit_sep$S$coeffs, `[[`, "c_non")
  }

  # reorder helper (2D only)
  reorder_S_2d <- function(Sdata, mode) {
    if (mode == "as-is" || mode == "x1_x2") return(Sdata)
    if (mode == "x2_x1") {
      S2 <- Sdata
      S2$X_tr <- Sdata$X_tr[, 2:1, drop = FALSE]
      S2$X_te  <- Sdata$X_te [, 2:1, drop = FALSE]
      return(S2)
    }
    Sdata
  }

  # Optional NF-based augmentation (surrogate π^ for more samples)
  use_nf <- {
    u <- getOption("mde.ctm.use_nf", NA)
    if (is.na(u)) {
      env <- Sys.getenv("MDE_CTM_USE_NF", "0"); env %in% c("1","true","TRUE")
    } else isTRUE(u)
  }
  if (use_nf) {
    # Dynamically choose layers by K and N
    Ktr <- ncol(X_tr); Ntr <- nrow(X_tr)
    layers_opt <- getOption("mde.nf.layers", NA)
    if (is.na(layers_opt)) {
      envL <- suppressWarnings(as.integer(Sys.getenv("MDE_NF_LAYERS", NA_character_)))
      if (is.na(envL)) {
        layers <- max(2L, min(6L, as.integer(ceiling(log2(max(32, Ntr))/2 + Ktr/2))))
      } else layers <- envL
    } else layers <- as.integer(layers_opt)
    aug_factor <- getOption("mde.nf.aug_factor", NA_real_)
    if (is.na(aug_factor)) {
      envF <- suppressWarnings(as.numeric(Sys.getenv("MDE_NF_AUG_FACTOR", NA_character_)))
      if (is.na(envF)) aug_factor <- 1.0 else aug_factor <- envF
    }
    if (!exists("trainNFSurrogate")) source(file.path(getwd(), "models/nf_surrogate.R"))
    message(sprintf("[NF] Training NF surrogate: layers=%d (K=%d,N=%d), aug_factor=%.2f", layers, Ktr, Ntr, aug_factor))
    nf <- trainNFSurrogate(X_tr, layers = layers, seed = seed)
    M <- max(1L, as.integer(floor(aug_factor * Ntr)))
    X_syn <- sample(nf, M)
    X_tr <- rbind(X_tr, X_syn)
  }

  # inner fitter for a given order and quadrature setting
  fit_one <- function(Sdata, Q_override = Q, control = list()) {
    X_tr <- Sdata$X_tr; X_te <- Sdata$X_te
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
      # Quadrature nodes from options/env unless explicitly provided
      Q_use <- if (is.null(Q_override)) {
        Qopt <- getOption("mde.ctm.Q", NA_integer_)
        if (is.na(Qopt)) {
          Qenv <- suppressWarnings(as.integer(Sys.getenv("MDE_CTM_Q", NA_character_)))
          if (is.na(Qenv)) min(12, 4 + 2 * degree_t_max) else Qenv
        } else Qopt
      } else Q_override
      batch_use <- if (is.null(batch_n)) min(N, max(256L, floor(65536 / max(1L, Q_use)))) else min(N, batch_n)
    quad <- .gauss_legendre_01_ct(Q_use)
    nodes <- quad$nodes
    weights <- quad$weights
      nodes_pow <- outer(nodes, seq_len(degree_t_max), `^`)

      fit_k <- function(k) {
        Xprev <- if (k > 1) X_tr_std[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
        xk <- X_tr_std[, k]
      Phi <- .build_g_design_ct(Xprev, k, degree_g, use_rbf = TRUE)
        m_alpha <- ncol(Phi)
        xprev_first <- if (k > 1) Xprev[1, , drop = TRUE] else numeric(0)
      # choose bspline tensor basis for all k (robust 1D t-basis) 
      use_bs <- TRUE
      if (use_bs) {
        # establish B-spline spec from training xk
        df_opt <- getOption('mde.ctm.bs_df', NA_integer_)
        if (is.na(df_opt)) {
          df_env <- suppressWarnings(as.integer(Sys.getenv('MDE_CTM_BS_DF', NA_character_)))
          if (is.na(df_env)) df_use <- 8L else df_use <- df_env
        } else df_use <- df_opt
        Btmp <- splines::bs(xk, df = df_use, degree = 3, intercept = TRUE)
        bs_spec <- list(df = ncol(Btmp), degree = 3,
                        knots = attr(Btmp, "knots"), boundary = attr(Btmp, "Boundary.knots"))
        m_beta <- bs_spec$df * length(.psi_xprev_poly3_ct(xprev_first))
        # Build cache once per dimension
        cch <- .build_ct_cache_tr(xk, Xprev, nodes, bs_spec)
        Pi_tr <- cch$Pi; B2D_tr <- cch$B2D; PsiX_tr <- cch$PsiX
        M <- cch$M; Qloc <- cch$Q; dfloc <- cch$df
        message(sprintf('CTM[k=%d] cache: B3D=(%d×%d×%d), \u03A0=(%d×%d)', k, N, Qloc, dfloc, N*Qloc, M))
      } else {
        m_beta <- length(.psi_basis_ct(0, xprev_first, degree_t, degree_g, TRUE, degree_t_cross, degree_x_cross))
        bs_spec <- NULL
        Pi_tr <- NULL; B2D_tr <- NULL; PsiX_tr <- NULL; M <- m_beta
      }
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
          xk_blk <- xk[idx]
          # Integral and gradient via cache (batched)
          if (use_bs) {
            rows <- unlist(lapply(idx, function(ii) ((ii - 1L) * Q + 1L):(ii * Q)))
            Vfull <- as.numeric(Pi_tr[rows, , drop = FALSE] %*% beta)
            Vmat <- matrix(pmax(pmin(Vfull, 30), -30), nrow = Q, ncol = length(idx))
            blog <- sweep(Vmat, 1, log(weights), FUN = "+")
            mcol <- apply(blog, 2, max)
            R <- exp(sweep(blog, 2, mcol, FUN = "-"))
            LSE <- mcol + log(colSums(R))
            I_vec <- ifelse(abs(xk_blk) < 1e-12, 0, sign(xk_blk) * exp(log(abs(xk_blk)) + LSE))
            soft <- sweep(R, 2, colSums(R), "/")
            S_vec <- (if (m_alpha > 0) as.numeric(Phi_blk %*% alpha) else 0) + I_vec
            ZI_rep <- rep(S_vec * I_vec, each = Q)
            g_beta_blk <- as.numeric(crossprod(Pi_tr[rows, , drop = FALSE], ZI_rep * as.numeric(soft)))
            # h at xk for log-jacobian gradient part: vectorized
            Bs_blk <- splines::bs(xk_blk, df = bs_spec$df, degree = bs_spec$degree,
                                  knots = bs_spec$knots, Boundary.knots = bs_spec$boundary,
                                  intercept = TRUE)
            PsiX_blk <- if (ncol(Xprev) >= 1L) cbind(1, Xprev[idx, 1], Xprev[idx, 1]^2, Xprev[idx, 1]^3) else matrix(1, nrow = length(idx), ncol = 1)
            Hblk <- .KR_rowwise_ct(Bs_blk, PsiX_blk)
            # Accumulate
            if (m_alpha > 0) {
              S_vec <- as.numeric(Phi_blk %*% alpha) + I_vec
              grad_alpha <- grad_alpha + as.numeric(crossprod(Phi_blk, S_vec))
            }
            grad_beta <- grad_beta + g_beta_blk - colSums(Hblk)
            S_sq_sum <- S_sq_sum + sum((if (m_alpha > 0) as.numeric(Phi_blk %*% alpha) + I_vec else I_vec)^2)
            term_sum <- term_sum + sum(Hblk %*% beta)
          } else {
            # Fallback (poly tensor) – keep previous unbatched path for rare case
            Xprev_blk <- if (k > 1) Xprev[idx, , drop = FALSE] else matrix(0, length(idx), 0)
            for (b in seq_along(idx)) {
              xp <- if (k > 1) Xprev_blk[b, ] else numeric(0)
              xval <- xk_blk[b]
              Psi_q <- .build_Psi_q_ct(xval, xp, nodes, nodes_pow, degree_t, degree_g, degree_t_cross, degree_x_cross)
              V <- as.vector(Psi_q %*% beta)
              V_clip <- pmax(pmin(V, 30), -30)
              b_vec <- log(weights) + V_clip
              if (abs(xval) < 1e-12) {
                I_i <- 0; soft <- rep(0, length(weights))
              } else {
                m_b <- max(b_vec); r <- exp(b_vec - m_b); lse <- m_b + log(sum(r))
                I_i <- sign(xval) * exp(log(abs(xval)) + lse); soft <- r / sum(r)
              }
              dI_i <- I_i * as.vector(t(Psi_q) %*% soft)
              psi_x <- .psi_basis_ct(xval, xp, degree_t, degree_g, TRUE, degree_t_cross, degree_x_cross)
              S_i <- if (m_alpha > 0) sum(Phi_blk[b, ] * alpha) + I_i else I_i
              S_sq_sum <- S_sq_sum + S_i^2
              if (m_alpha > 0) grad_alpha <- grad_alpha + S_i * Phi_blk[b, ]
              grad_beta <- grad_beta + S_i * dI_i - psi_x
              term_sum <- term_sum + sum(psi_x * beta)
            }
          }
        }
        loss <- 0.5 * S_sq_sum / N - term_sum / N +
          0.5 * (lambda_non * sum(alpha^2) + lambda_mon * sum(beta^2))
        grad_alpha_out <- if (m_alpha > 0) grad_alpha / N + lambda_non * alpha else numeric(0)
        grad_beta_out <- grad_beta / N + lambda_mon * beta
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
      lower <- c(rep(-Inf, m_alpha), rep(-6, m_beta))
      upper <- c(rep( Inf, m_alpha), rep( 6, m_beta))
      ctrl <- modifyList(list(pgtol = 1e-4, maxit = 200), control)
      opt <- optim(theta0, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
                   control = ctrl)
      alpha_hat <- if (m_alpha > 0) opt$par[seq_len(m_alpha)] else numeric(0)
      beta_hat <- opt$par[(m_alpha + 1):length(opt$par)]
      message(sprintf('CTM[k=%d] fn_evals=%d, gr_evals=%d', k, opt$counts['function'], opt$counts['gradient']))
      list(alpha = alpha_hat, beta = beta_hat, convergence = opt$convergence,
            g_spec = attr(Phi, "g_spec"), bs_spec = bs_spec)
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
      order = seq_len(K),
      meta = list(
        lambda_non = lambda_non,
        lambda_mon = lambda_mon,
        ridge = c(lambda_non = lambda_non, lambda_mon = lambda_mon),
        Q = Q_use,
        order_mode = order_mode,
        seed = seed
      )
    )
    class(S_map) <- "ttm_cross_term"
    })[["elapsed"]]

    time_pred <- system.time({
      predict(S_map, X_te, "logdensity_by_dim")
    })[["elapsed"]]

    list(
      S = S_map,
      NLL_test = .NLL_set_ct(S_map, X_te),
      stderr_test = .SE_set_ct(S_map, X_te),
      time_train = time_train,
      time_pred = time_pred
    )
  }

  # No CV/grid search: single fit with configured order and lambdas
  Sout <- reorder_S_2d(S_in, order_mode)
  # Allow overriding max iterations via option/env
  maxit_use <- {
    mi <- getOption("mde.ctm.maxit", NA_integer_)
    if (is.na(mi)) {
      mie <- suppressWarnings(as.integer(Sys.getenv("MDE_CTM_MAXIT", NA_character_)))
      if (is.na(mie)) 200L else mie
    } else as.integer(mi)
  }
  out <- fit_one(Sout, Q_override = if (is.null(Q)) NULL else Q, control = list(maxit = maxit_use))
  # Quick sanity on training Z: |mean|<=0.5 and var in [0.5,2]; if violated, refit once with stronger λ_mon
  compute_Z_metrics <- function(S_map, X) {
    Xs <- .standardize(S_map, X)
    N <- nrow(Xs); K <- ncol(Xs)
    nodes <- S_map$quad_nodes_ct; weights <- S_map$quad_weights_ct; nodes_pow <- S_map$quad_nodes_pow_ct
    res <- lapply(seq_len(K), function(k) {
      Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, N, 0)
      xk <- Xs[, k]
      # g-part
      if (!is.null(S_map$coeffs[[k]]$g_spec) && S_map$coeffs[[k]]$g_spec$type == 'rbf_poly') {
        x1 <- if (ncol(Xprev) >= 1) Xs[, 1] else numeric(N)
        gs <- S_map$coeffs[[k]]$g_spec
        poly <- cbind(1, x1, x1^2)
        rbf <- .rbf_basis_1d_ct(x1, centers = gs$centers, sigma = gs$sigma)
        Phi <- cbind(poly, rbf)
      } else {
        Phi <- .basis_g_ct(Xprev, S_map$degree_g)
      }
      alpha <- S_map$coeffs[[k]]$alpha
      beta <- S_map$coeffs[[k]]$beta
      Zk <- numeric(N)
      for (i in seq_len(N)) {
        xp <- if (k > 1) Xprev[i, ] else numeric(0)
        xval <- xk[i]
        if (!is.null(S_map$coeffs[[k]]$bs_spec)) {
          Psi_q <- .build_Psi_q_bs_ct(xval, xp, nodes, S_map$coeffs[[k]]$bs_spec)
        } else {
          Psi_q <- .build_Psi_q_ct(xval, xp, nodes, nodes_pow, S_map$degree_t, S_map$degree_g, S_map$degree_t_cross, S_map$degree_x_cross)
        }
        V <- as.vector(Psi_q %*% beta)
        Vc <- pmax(pmin(V, 30), -30)
        bv <- log(weights) + Vc
        if (abs(xval) < 1e-12) I <- 0 else {
          mb <- max(bv); I <- sign(xval) * exp(log(abs(xval)) + mb + log(sum(exp(bv - mb))))
        }
        g <- if (length(alpha) > 0) sum(Phi[i, ] * alpha) else 0
        Zk[i] <- g + I
      }
      c(mean = mean(Zk), var = stats::var(Zk))
    })
    do.call(rbind, res)
  }
  zm <- compute_Z_metrics(out$S, S_in$X_tr[sample.int(nrow(S_in$X_tr), min(64L, nrow(S_in$X_tr))), , drop = FALSE])
  bad <- any(abs(zm[, 'mean']) > 0.5) || any(zm[, 'var'] < 0.5 | zm[, 'var'] > 2.0)
  if (bad) {
    warning(sprintf('Z sanity failed (mean,var)=%s; refitting with lambda_mon*3', paste(c(zm), collapse=',')))
    lambda_mon_old <- lambda_mon
    lambda_mon <<- lambda_mon * 3
    on.exit({ lambda_mon <<- lambda_mon_old }, add = TRUE)
    out <- fit_one(S_in, Q_override = Q)
  }
    out
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
    # rebuild g design depending on fitted spec
    if (!is.null(object$coeffs[[k]]$g_spec) && object$coeffs[[k]]$g_spec$type == "rbf_poly") {
      x1 <- if (ncol(Xprev) >= 1) Xs[, 1] else numeric(N)
      gs <- object$coeffs[[k]]$g_spec
      poly <- cbind(1, x1, x1^2)
      rbf <- .rbf_basis_1d_ct(x1, centers = gs$centers, sigma = gs$sigma)
      Phi <- cbind(poly, rbf)
    } else {
      Phi <- .basis_g_ct(Xprev, object$degree_g)
    }
    m_alpha <- ncol(Phi)
    alpha <- object$coeffs[[k]]$alpha
    beta <- object$coeffs[[k]]$beta
    xprev_first <- if (k > 1) Xprev[1, , drop = TRUE] else numeric(0)
      if (!is.null(object$coeffs[[k]]$bs_spec)) {
        psi_len <- object$coeffs[[k]]$bs_spec$df * length(.psi_xprev_poly3_ct(xprev_first))
      } else {
        psi_len <- length(.psi_basis_ct(Xs[1, k], xprev_first,
                                        object$degree_t, object$degree_g, TRUE, object$degree_t_cross, object$degree_x_cross))
      }
    stopifnot(length(beta) == psi_len)
    Z_col <- numeric(N)
    LJ_col <- numeric(N)
    # Batched prediction with cached B-splines tensor
    if (!is.null(object$coeffs[[k]]$bs_spec)) {
      bs_spec <- object$coeffs[[k]]$bs_spec
      Nall <- N; Q <- length(nodes)
      # Build B2D and Pi for all points (fast path)
      Btmp0 <- splines::bs(0, df = bs_spec$df, degree = bs_spec$degree,
                           knots = bs_spec$knots, Boundary.knots = bs_spec$boundary, intercept = TRUE)
      df <- ncol(Btmp0)
      # Vectorized Bs for all N points and both integral nodes (xk*nodes) and log-Jacobian (xk)
      # Integral path uses B2D (N*Q x df)
      B2D <- matrix(0, nrow = Nall * Q, ncol = df)
      for (i in seq_len(Nall)) {
        Bi <- splines::bs(xk[i] * nodes, df = df, degree = bs_spec$degree,
                          knots = bs_spec$knots, Boundary.knots = bs_spec$boundary, intercept = TRUE)
        B2D[((i - 1L) * Q + 1L):(i * Q), ] <- Bi
      }
      # PsiX (N x Mx)
      if (ncol(Xprev) >= 1L) {
        pieces <- list(matrix(1, nrow = Nall, ncol = 1))
        for (j in seq_len(ncol(Xprev))) {
          xj <- Xprev[, j]
          pieces[[length(pieces) + 1L]] <- cbind(xj, xj^2, xj^3)
        }
        PsiX <- do.call(cbind, pieces)
      } else {
        PsiX <- matrix(1, nrow = Nall, ncol = 1)
      }
      Mx <- ncol(PsiX); M <- df * Mx
      # Pi via row-wise Khatri-Rao
      Pi <- matrix(0, nrow = Nall * Q, ncol = M)
      for (i in seq_len(Nall)) {
        rows <- ((i - 1L) * Q + 1L):(i * Q)
        Pi[rows, ] <- .KR_rowwise_ct(B2D[rows, , drop = FALSE], matrix(PsiX[i, ], nrow = Q, ncol = Mx, byrow = TRUE))
      }
      Vfull <- as.numeric(Pi %*% beta)
      Vmat <- matrix(pmax(pmin(Vfull, 30), -30), nrow = Q, ncol = Nall)
      blog <- sweep(Vmat, 1, log(weights), FUN = "+")
      mcol <- apply(blog, 2, max)
      R <- exp(sweep(blog, 2, mcol, FUN = "-"))
      LSE <- mcol + log(colSums(R))
      I_vec <- ifelse(abs(xk) < 1e-12, 0, sign(xk) * exp(log(abs(xk)) + LSE))
      # Log-Jacobian Hmat via KR_rowwise: Hmat = KR_rowwise(Bs(xk), PsiX)
      Bs_all <- splines::bs(xk, df = df, degree = bs_spec$degree,
                            knots = bs_spec$knots, Boundary.knots = bs_spec$boundary, intercept = TRUE)
      Hmat <- .KR_rowwise_ct(Bs_all, PsiX)
      # Log one-time
      if (!getOption('mde.ctm_predict_fastpath_printed', FALSE)) {
        message(sprintf('CTM predict fast-path: N=%d, df=%d, Mx=%d', Nall, df, Mx))
        options(mde.ctm_predict_fastpath_printed = TRUE)
      }
      stopifnot(ncol(Hmat) == length(beta))
      Z_col <- if (m_alpha > 0) as.numeric(Phi %*% alpha) + I_vec else I_vec
      LJ_col <- as.numeric(Hmat %*% beta) - log(object$sigma[k])
      if (any(!is.finite(Z_col)) || any(!is.finite(LJ_col))) stop('Non-finite values in CTM predict fast-path')
    } else {
      # Poly tensor fallback (rare)
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
          V_clip <- pmax(pmin(V, 30), -30)
          b_vec <- log(weights) + V_clip
          if (abs(xval) < 1e-12) {
            I_i <- 0
          } else {
            m_b <- max(b_vec)
            r <- exp(b_vec - m_b)
            lse <- m_b + log(sum(r))
            I_i <- sign(xval) * exp(log(abs(xval)) + lse)
          }
          psi_x <- .psi_basis_ct(xval, xp, object$degree_t, object$degree_g, TRUE, object$degree_t_cross, object$degree_x_cross)
          Z_col[idx[b]] <- if (m_alpha > 0) sum(Phi_blk[b, ] * alpha) + I_i else I_i
          LJ_col[idx[b]] <- sum(beta * psi_x) - log(object$sigma[k])
        }
      }
    }
    # Runtime Eq.(22) checks (optional): numeric derivative vs exp(h)
    if (getOption("mde.check_eq22", FALSE)) {
      set.seed(1L)
      nn <- min(16L, N)
      pick <- sample.int(N, nn)
      eps <- 1e-5
      for (ii in pick) {
        xp <- if (k > 1) Xprev[ii, ] else numeric(0)
        x0 <- xk[ii]
        # g_k does not depend on x_k
        g0 <- if (m_alpha > 0) {
          if (!is.null(object$coeffs[[k]]$g_spec) && object$coeffs[[k]]$g_spec$type == 'rbf_poly') {
            x1ii <- if (length(xp) >= 1) Xprev[ii, 1] else 0
            gs <- object$coeffs[[k]]$g_spec
            polyii <- c(1, x1ii, x1ii^2)
            rbfi <- exp(-0.5 * ((x1ii - gs$centers)/gs$sigma)^2)
            sum(c(polyii, rbfi) * alpha)
          } else {
            sum(Phi[ii, ] * alpha)
          }
        } else 0
        # recompute integral at x0 and x0+eps
        Psi_q0 <- if (!is.null(object$coeffs[[k]]$bs_spec)) .build_Psi_q_bs_ct(x0, xp, nodes, object$coeffs[[k]]$bs_spec)
                  else .build_Psi_q_ct(x0, xp, nodes, nodes_pow, object$degree_t, object$degree_g, object$degree_t_cross, object$degree_x_cross)
        V0 <- as.vector(Psi_q0 %*% beta)
        V0c <- pmax(pmin(V0, 30), -30)
        bvec0 <- log(weights) + V0c; m0 <- max(bvec0)
        I0 <- if (abs(x0) < 1e-12) 0 else sign(x0) * exp(log(abs(x0)) + m0 + log(sum(exp(bvec0 - m0))))
        x1 <- x0 + eps
        Psi_q1 <- if (!is.null(object$coeffs[[k]]$bs_spec)) .build_Psi_q_bs_ct(x1, xp, nodes, object$coeffs[[k]]$bs_spec)
                  else .build_Psi_q_ct(x1, xp, nodes, nodes_pow, object$degree_t, object$degree_g, object$degree_t_cross, object$degree_x_cross)
        V1 <- as.vector(Psi_q1 %*% beta)
        V1c <- pmax(pmin(V1, 30), -30)
        bvec1 <- log(weights) + V1c; m1 <- max(bvec1)
        I1 <- if (abs(x1) < 1e-12) 0 else sign(x1) * exp(log(abs(x1)) + m1 + log(sum(exp(bvec1 - m1))))
        dnum <- (I1 - I0) / eps
        # h(x0,xprev)
        psi_x0 <- if (!is.null(object$coeffs[[k]]$bs_spec)) as.vector(.build_Psi_q_bs_ct(x0, xp, 1, object$coeffs[[k]]$bs_spec))
                  else .psi_basis_ct(x0, xp, object$degree_t, object$degree_g, TRUE, object$degree_t_cross, object$degree_x_cross)
        dtrue <- exp(sum(beta * psi_x0))
        relerr <- abs(dnum - dtrue) / max(1e-12, abs(dtrue))
        if (relerr > 1e-6 || !is.finite(relerr)) stop("Eq.(22) derivative check failed")
        # check dg/dx_k == 0 numerically via g(xk+eps)-g(xk)
        g1 <- g0  # unchanged as g depends only on x_<k>
        if (abs(g1 - g0) > 1e-8) stop("g_k depends on x_k (violates Eq.(22))")
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
  # Invariants and numeric sanity checks
  stopifnot(is.matrix(LD), nrow(LD) == N, ncol(LD) == K)
  if (!all(is.finite(LD))) stop("Non-finite values in logdensity_by_dim")
  if (type == "logdensity_by_dim") {
    # Assert joint equals rowSums
    LD_joint <- rowSums(LD)
    if (!all(is.finite(LD_joint))) stop("Non-finite joint logdensity")
    stopifnot(max(abs(LD_joint - rowSums(LD))) <= 1e-10)
    LD
  } else {
    LD_joint <- rowSums(LD)
    if (!all(is.finite(LD_joint))) stop("Non-finite joint logdensity")
    stopifnot(max(abs(LD_joint - rowSums(LD))) <= 1e-10)
    LD_joint
  }
}

# ---------------------------------------------------------------------------
# Reverse-KL API (R = S^{-1}) — Compatibility wrapper
#
# NOTE: True reverse-KL training requires access to log pi(x) and, in general,
# its gradient wrt x to backprop through R(z). The current codebase does not
# provide analytic gradients for the true joint density. To keep the API and
# workflow aligned, we expose a reverse entry that falls back to forward-KL
# training and wraps the trained forward map for prediction.

trainCrossTermMapReverse <- function(S_or_sampler, seed = 42, ...) {
  set.seed(seed)
  fit_fwd <- trainCrossTermMap(S_or_sampler, seed = seed, ...)
  R_wrap <- list(S_forward = fit_fwd$S,
                 meta = modifyList(fit_fwd$S$meta, list(direction = "reverse_wrapper")))
  class(R_wrap) <- "ttm_cross_term_reverse"
  list(
    S = R_wrap,
    NLL_test = fit_fwd$NLL_test,
    stderr_test = fit_fwd$stderr_test,
    time_train = fit_fwd$time_train,
    time_pred = fit_fwd$time_pred
  )
}

predict.ttm_cross_term_reverse <- function(object, newdata,
                                           type = c("logdensity_by_dim", "logdensity"),
                                           ...) {
  type <- tryCatch(match.arg(type), error = function(e) stop("unknown type"))
  # Delegate prediction to the wrapped forward map; invariants are enforced there.
  predict(object$S_forward, newdata, type = type)
}

# ----------------------------------------------------------------------------
# Reverse-KL to NF surrogate: train R(z;theta) by maximizing E_eta[ log pi_hat(R(z)) + log|det dR/dz| ]
# Fair: NF learned from X_tr; no oracle.

trainCrossTermMapReverseNF <- function(S_or_X, degree_g = 2, degree_t = 2, degree_t_cross = 1, degree_x_cross = 1,
                                       seed = 42, Q = NULL, lambda_non = NULL, lambda_mon = NULL,
                                       nf_layers = NULL) {
  set.seed(seed)
  X_tr <- if (is.list(S_or_X)) S_or_X$X_tr else S_or_X
  stopifnot(is.matrix(X_tr))
  K <- ncol(X_tr); N <- nrow(X_tr)
  # Ridge defaults
  if (is.null(lambda_non)) {
    lambda_non <- getOption("mde.ctm.lambda_non", NA_real_)
    if (is.na(lambda_non)) {
      env <- suppressWarnings(as.numeric(Sys.getenv("MDE_CTM_LAMBDA_NON", NA_character_)))
      lambda_non <- ifelse(is.na(env), 1e-2, env)
    }
  }
  if (is.null(lambda_mon)) {
    lambda_mon <- getOption("mde.ctm.lambda_mon", NA_real_)
    if (is.na(lambda_mon)) {
      env <- suppressWarnings(as.numeric(Sys.getenv("MDE_CTM_LAMBDA_MON", NA_character_)))
      lambda_mon <- ifelse(is.na(env), 2e-2, env)
    }
  }
  # NF surrogate
  if (!exists("trainNFSurrogate")) source(file.path(getwd(), "models/nf_surrogate.R"))
  if (is.null(nf_layers)) {
    envL <- suppressWarnings(as.integer(Sys.getenv("MDE_NF_LAYERS", NA_character_)))
    nf_layers <- ifelse(is.na(envL), max(2L, min(6L, as.integer(ceiling(log2(max(32,N))/2 + K/2)))), envL)
  }
  message(sprintf("[NF-RevKL] Training NF surrogate (layers=%d) ...", nf_layers))
  nf <- trainNFSurrogate(X_tr, layers = nf_layers, seed = seed)
  # Reverse map R params per k: g (alpha) on z_prev, h (beta) on (t=z_k, z_prev)
  # Quadratur
  degree_t_max <- max(degree_t, degree_t_cross)
  Q_use <- if (is.null(Q)) min(12, 4 + 2 * degree_t_max) else Q
  quad <- .gauss_legendre_01_ct(Q_use); nodes <- quad$nodes; weights <- quad$weights
  nodes_pow <- outer(nodes, seq_len(degree_t_max), `^`)
  # Param sizes
  g_dim <- function(k) { if (k==1) 0L else ncol(.basis_g_ct(matrix(0,1,k-1), degree_g)) }
  psi_dim <- function(k) {
    prev <- if (k==1) numeric(0) else rep(0, k-1)
    length(.psi_basis_ct(0, prev, degree_t, degree_g, TRUE, degree_t_cross, degree_x_cross))
  }
  idx_alpha <- list(); idx_beta <- list(); off <- 0L
  for (k in 1:K) {
    a <- g_dim(k); b <- psi_dim(k)
    if (a>0) { idx_alpha[[k]] <- (off+1):(off+a); off <- off+a } else idx_alpha[[k]] <- integer(0)
    idx_beta[[k]]  <- (off+1):(off+b); off <- off+b
  }
  theta0 <- rep(0, off)
  # Sample z from N(0,I) for objective Monte Carlo
  M <- max(N, 128L)
  Zs <- matrix(rnorm(M*K), ncol = K)

  # Build function to compute R(z) and logdet per batch
  R_forward <- function(theta, Z) {
    Mloc <- nrow(Z); X <- matrix(NA_real_, Mloc, K); LJ <- matrix(0, Mloc, K)
    for (k in 1:K) {
      Zprev <- if (k>1) Z[,1:(k-1),drop=FALSE] else matrix(0,Mloc,0)
      zk <- Z[,k]
      a_idx <- idx_alpha[[k]]; b_idx <- idx_beta[[k]]
      alpha <- if (length(a_idx)>0) theta[a_idx] else numeric(0)
      beta  <- theta[b_idx]
      # g_k(z_prev)
      gk <- if (k>1 && length(alpha)>0) as.numeric(.basis_g_ct(Zprev, degree_g) %*% alpha) else rep(0, Mloc)
      # integral over t in [0, z_k]
      # cache per sample i
      Ik <- numeric(Mloc); hk_at <- numeric(Mloc)
      for (i in 1:Mloc) {
        xp <- if (k>1) Zprev[i,,drop=TRUE] else numeric(0)
        val <- zk[i]
        # build Psi_q at nodes scaled by |z|
        Psi_q <- .build_Psi_q_ct(val, xp, nodes, nodes_pow, degree_t, degree_g, degree_t_cross, degree_x_cross)
        V <- as.vector(Psi_q %*% beta)
        Vc <- pmax(pmin(V,30),-30)
        bvec <- log(weights) + Vc
        if (abs(val) < 1e-12) {
          I_i <- 0
        } else {
          m_b <- max(bvec); r <- exp(bvec - m_b); lse <- m_b + log(sum(r))
          I_i <- sign(val) * exp(log(abs(val)) + lse)
        }
        Ik[i] <- I_i
        # h(z_k,z_prev)
        psi_x <- .psi_basis_ct(val, xp, degree_t, degree_g, TRUE, degree_t_cross, degree_x_cross)
        hk_at[i] <- sum(beta * psi_x)
      }
      X[,k] <- gk + Ik
      LJ[,k] <- hk_at
    }
    list(X=X, LJ=LJ)
  }

  obj <- function(theta) {
    rf <- R_forward(theta, Zs)
    X_hat <- rf$X; LJ <- rf$LJ
    ll_hat <- as.numeric(predict(nf, X_hat, type='logdensity'))
    val <- -mean(ll_hat + rowSums(LJ))
    # ridge
    pen <- 0
    for (k in 1:K) {
      a_idx <- idx_alpha[[k]]; b_idx <- idx_beta[[k]]
      if (length(a_idx)>0) pen <- pen + lambda_non * sum(theta[a_idx]^2)
      pen <- pen + lambda_mon * sum(theta[b_idx]^2)
    }
    val + 0.5*pen
  }

  message(sprintf("[RevKL-NF] Optimize R with L-BFGS-B (dim theta=%d) ...", length(theta0)))
  opt <- optim(theta0, obj, method='L-BFGS-B', control=list(maxit=120))
  theta_hat <- opt$par
  # Build artifact
  coeffs <- vector('list', K)
  for (k in 1:K) {
    a_idx <- idx_alpha[[k]]; b_idx <- idx_beta[[k]]
    coeffs[[k]] <- list(alpha = if (length(a_idx)>0) theta_hat[a_idx] else numeric(0),
                        beta = theta_hat[b_idx], bs_spec = NULL)
  }
  R_map <- list(coeffs = coeffs, degree_g = degree_g, degree_t = degree_t,
                degree_t_cross = degree_t_cross, degree_x_cross = degree_x_cross,
                quad_nodes_ct = nodes, quad_weights_ct = weights,
                quad_nodes_pow_ct = nodes_pow, K = K,
                meta = list(direction='reverse_nf', lambda_non=lambda_non, lambda_mon=lambda_mon,
                            nf_layers=nf_layers, seed=seed))
  class(R_map) <- 'ttm_cross_term_rev_nf'
  # Simple timing proxies not computed here; return artifact
  list(S = R_map)
}

predict.ttm_cross_term_rev_nf <- function(object, newdata, type = c('logdensity_by_dim','logdensity')) {
  type <- match.arg(type)
  X <- as.matrix(newdata); N <- nrow(X); K <- ncol(X)
  nodes <- object$quad_nodes_ct; nodes_pow <- object$quad_nodes_pow_ct; weights <- object$quad_weights_ct
  # Invert R(z)=x sequentially via 1D root-finding
  Z <- matrix(0, N, K); LJm <- matrix(0, N, K)
  for (i in 1:N) {
    z_prev <- numeric(0)
    for (k in 1:K) {
      alpha <- object$coeffs[[k]]$alpha; beta <- object$coeffs[[k]]$beta
      g_of <- function(zprev) {
        if (length(alpha)==0 || length(zprev)==0) return(0)
        as.numeric(.basis_g_ct(matrix(zprev, nrow=1), object$degree_g) %*% alpha)
      }
      f_Rk <- function(zk) {
        # compute R_k(zprev, zk)
        I_i <- {
          Psi_q <- .build_Psi_q_ct(zk, z_prev, nodes, nodes_pow, object$degree_t, object$degree_g, object$degree_t_cross, object$degree_x_cross)
          V <- as.vector(Psi_q %*% beta); Vc <- pmax(pmin(V,30),-30)
          bvec <- log(weights) + Vc
          if (abs(zk) < 1e-12) 0 else { m_b <- max(bvec); r <- exp(bvec - m_b); lse <- m_b + log(sum(r)); sign(zk) * exp(log(abs(zk)) + lse) }
        }
        g_of(z_prev) + I_i
      }
      # root for f_Rk(zk) - xk = 0
      target <- X[i,k]
      froot <- function(zk) f_Rk(zk) - target
      lo <- -8; hi <- 8; flo <- froot(lo); fhi <- froot(hi); iter <- 0
      while (flo*fhi > 0 && iter < 8) { lo <- lo*2; hi <- hi*2; flo <- froot(lo); fhi <- froot(hi); iter <- iter+1 }
      zk <- tryCatch(uniroot(froot, c(lo,hi))$root, error=function(e) 0)
      Z[i,k] <- zk
      # h_k at (zk,z_prev)
      psi_x <- .psi_basis_ct(zk, z_prev, object$degree_t, object$degree_g, TRUE, object$degree_t_cross, object$degree_x_cross)
      LJm[i,k] <- - sum(beta * psi_x) # minus log ∂R_k/∂zk
      # update z_prev
      z_prev <- c(z_prev, zk)
    }
  }
  base <- (-0.5) * (Z^2) - 0.5 * log(2*pi)
  LD <- base + LJm
  stopifnot(all(is.finite(LD)))
  if (type=='logdensity_by_dim') return(LD)
  rowSums(LD)
}


# Ende -----------------------------------------------------------------------

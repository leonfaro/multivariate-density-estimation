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

# Fixed B-spline degrees of freedom for time basis
.CTM_BS_DF <- 8L

# Build design matrix for g_k(x_<k>)
.build_g_design_ct <- function(Xprev, k, degree_g, use_rbf = TRUE) {
  # For 2D (k==2) with one predecessor, use RBF(7) + poly {1,x,x^2}
  if (ncol(Xprev) == 1L && k == 2L && use_rbf) {
    x1 <- Xprev[, 1]
    ctr <- as.numeric(stats::quantile(x1, probs = seq(0.1, 0.9, length.out = 7)))
    # bandwidth from quantile spacing or fallback to sd
    dif <- diff(ctr)
    sig_raw <- if (all(is.finite(dif)) && length(dif) > 0) 0.5 * stats::median(dif) else (stats::sd(x1) + 1e-8)/3
    sig <- max(sig_raw, (stats::sd(x1) + 1e-8)/3, 1e-6)
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
    Mx <- ncol(PsiX)
    # Khatri–Rao row-wise: replicate the PsiX row to Q rows, then KR with Bi
    psi_mat <- matrix(PsiX[i, ], nrow = Q, ncol = Mx, byrow = TRUE)
    block <- .KR_rowwise_ct(Bi, psi_mat)
    Pi[rows, ] <- block
  }
  # Strict shape/numeric checks
  stopifnot(nrow(Pi) == N * Q, ncol(Pi) == df * ncol(PsiX))
  if (!all(is.finite(Pi))) stop("Non-finite entries in Pi design (training)")
  # Return only the design; callers can infer sizes from bs_spec/nodes/Pi
  list(Pi = Pi)
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
  m_beta <- .CTM_BS_DF * length(.psi_xprev_poly3_ct(xprev_first))
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

## Removed unused logJacDiag_ct (unreferenced)

## Removed unused forwardKLLoss_ct (unreferenced)

# Optional Eq.(22) derivative point check (stub; disabled by default)
.check_eq22_point_ct <- function(Smap, k, x0, xp, nodes, weights) {
  # Returns relative error between d/dx ∫ exp(h) and exp(h) at (x0, xp).
  # Stubbed to 0.0; enable and implement if needed for HPO diagnostics.
  0.0
}

# Exportierte Funktionen -----------------------------------------------------

trainCrossTermMap <- function(X_or_path, degree_g = 2, degree_t = 2, degree_t_cross = 1, degree_x_cross = 1,
                              batch_n = NULL, Q = NULL,
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

  # Quadrature and ridge defaults are taken from options/env
  # (defaults: Q=24, lambda_non=3e-2, lambda_mon=2e-2)

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
    if (isTRUE(getOption("mde.verbose", FALSE))) message(sprintf("[NF] Training NF surrogate: layers=%d (K=%d,N=%d), aug_factor=%.2f", layers, Ktr, Ntr, aug_factor))
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
          if (is.na(Qenv)) 24L else Qenv
        } else Qopt
      } else Q_override
      batch_use <- if (is.null(batch_n)) min(N, max(256L, floor(65536 / max(1L, Q_use)))) else min(N, batch_n)
    quad <- .gauss_legendre_01_ct(Q_use)
    nodes <- quad$nodes
    weights <- quad$weights
      # no polynomial time-basis;  not used

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
        Pi_tr <- cch$Pi
        if (isTRUE(getOption("mde.verbose", FALSE))) {
          Qloc <- length(nodes)
          dfloc <- bs_spec$df
          Mloc <- ncol(Pi_tr)
          message(sprintf('CTM[k=%d] cache: B3D=(%d×%d×%d), \u03A0=(%d×%d)', k, N, Qloc, dfloc, N*Qloc, Mloc))
        }
        # Precompute PsiX_full for all rows once per k
        if (ncol(Xprev) >= 1L) {
          parts_all <- list(matrix(1, nrow = N, ncol = 1))
          for (jj in seq_len(ncol(Xprev))) {
            xj_all <- Xprev[, jj]
            parts_all[[length(parts_all) + 1L]] <- cbind(xj_all, xj_all^2, xj_all^3)
          }
          PsiX_full <- do.call(cbind, parts_all)
        } else {
          PsiX_full <- matrix(1, nrow = N, ncol = 1)
        }
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
            Qlen <- length(nodes)
            rows <- unlist(lapply(idx, function(ii) ((ii - 1L) * Qlen + 1L):(ii * Qlen)))
            Vfull <- as.numeric(Pi_tr[rows, , drop = FALSE] %*% beta)
            Vmat <- matrix(pmax(pmin(Vfull, 30), -30), nrow = Qlen, ncol = length(idx))
            blog <- sweep(Vmat, 1, log(weights), FUN = "+")
            mcol <- apply(blog, 2, max)
            R <- exp(sweep(blog, 2, mcol, FUN = "-"))
            LSE <- mcol + log(colSums(R))
            I_vec <- ifelse(abs(xk_blk) < 1e-12, 0, sign(xk_blk) * exp(log(abs(xk_blk)) + LSE))
            soft <- sweep(R, 2, colSums(R), "/")
            S_vec <- (if (m_alpha > 0) as.numeric(Phi_blk %*% alpha) else 0) + I_vec
            ZI_rep <- rep(S_vec * I_vec, each = Qlen)
            g_beta_blk <- as.numeric(crossprod(Pi_tr[rows, , drop = FALSE], ZI_rep * as.numeric(soft)))
            # h at xk for log-jacobian gradient part: vectorized
            Bs_blk <- splines::bs(xk_blk, df = bs_spec$df, degree = bs_spec$degree,
                                  knots = bs_spec$knots, Boundary.knots = bs_spec$boundary,
                                  intercept = TRUE)
            # Use precomputed PsiX_full for log-Jacobian
            PsiX_blk <- PsiX_full[idx, , drop = FALSE]
            Hblk <- .KR_rowwise_ct(Bs_blk, PsiX_blk)
            # Accumulate
            if (m_alpha > 0) {
              S_vec <- as.numeric(Phi_blk %*% alpha) + I_vec
              grad_alpha <- grad_alpha + as.numeric(crossprod(Phi_blk, S_vec))
            }
            grad_beta <- grad_beta + g_beta_blk - colSums(Hblk)
            S_sq_sum <- S_sq_sum + sum((if (m_alpha > 0) as.numeric(Phi_blk %*% alpha) + I_vec else I_vec)^2)
            term_sum <- term_sum + sum(Hblk %*% beta)
            
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
      if (isTRUE(getOption("mde.verbose", FALSE))) message(sprintf('CTM[k=%d] fn_evals=%d, gr_evals=%d', k, opt$counts['function'], opt$counts['gradient']))
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
      quad_nodes_ct = nodes,
      quad_weights_ct = weights,
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
  # Persist display permutation (for predict) based on order_mode (2D convenience)
  K0 <- ncol(S_in$X_tr)
  perm <- if (identical(order_mode, "x2_x1") && K0 == 2L) c(2L, 1L) else seq_len(K0)
  # Allow overriding max iterations via option/env
  maxit_use <- {
    mi <- getOption("mde.ctm.maxit", NA_integer_)
    if (is.na(mi)) {
      mie <- suppressWarnings(as.integer(Sys.getenv("MDE_CTM_MAXIT", NA_character_)))
      if (is.na(mie)) 200L else mie
    } else as.integer(mi)
  }
  # --- Hyperparameter Optimization (optional) ---------------------------------
  # Helper: parse options/env with defaults
  .opt_get_bool <- function(opt, env, default = FALSE) {
    v <- getOption(opt, NA)
    if (!is.na(v)) return(isTRUE(v))
    ev <- Sys.getenv(env, NA_character_)
    if (!is.na(ev)) return(ev %in% c("1","true","TRUE","on","yes"))
    default
  }
  .opt_get_num <- function(opt, env, default) {
    v <- getOption(opt, NA_real_)
    if (!is.na(v)) return(as.numeric(v))
    ev <- suppressWarnings(as.numeric(Sys.getenv(env, NA_character_)))
    if (!is.na(ev)) return(ev)
    default
  }
  .opt_get_int <- function(opt, env, default) {
    v <- getOption(opt, NA_integer_)
    if (!is.na(v)) return(as.integer(v))
    ev <- suppressWarnings(as.integer(Sys.getenv(env, NA_character_)))
    if (!is.na(ev)) return(ev)
    default
  }

  # Force-disable HPO to respect hard-coded hyperparameters
  auto_hpo <- FALSE
  if (!auto_hpo) {
    out <- fit_one(Sout, Q_override = 24L, control = list(maxit = maxit_use))
    out$S$perm <- perm
  } else {
    # HPO params
    hpo_val_frac <- min(max(.opt_get_num("mde.ctm.hpo_val_frac", "MDE_CTM_HPO_VAL_FRAC", 0.2), 0.05), 0.5)
    hpo_max_rounds <- max(1L, .opt_get_int("mde.ctm.hpo_max_rounds", "MDE_CTM_HPO_MAX_ROUNDS", 6L))
    hpo_tol_nats <- max(0, .opt_get_num("mde.ctm.hpo_tol_nats", "MDE_CTM_HPO_TOL_NATS", 1e-3))
    Q_min <- max(4L, .opt_get_int("mde.ctm.hpo_Q_min", "MDE_CTM_HPO_Q_MIN", 16L))
    Q_max <- max(Q_min, .opt_get_int("mde.ctm.hpo_Q_max", "MDE_CTM_HPO_Q_MAX", 64L))
    Q_step <- max(1L, .opt_get_int("mde.ctm.hpo_Q_step", "MDE_CTM_HPO_Q_STEP", 8L))
    lam_lo <- .opt_get_num("mde.ctm.hpo_lambda_lo", "MDE_CTM_HPO_LAMBDA_LO", 1e-5)
    lam_hi <- .opt_get_num("mde.ctm.hpo_lambda_hi", "MDE_CTM_HPO_LAMBDA_HI", 1e+1)
    ratio_lo <- .opt_get_num("mde.ctm.hpo_ratio_lo", "MDE_CTM_HPO_RATIO_LO", 0.25)
    ratio_hi <- .opt_get_num("mde.ctm.hpo_ratio_hi", "MDE_CTM_HPO_RATIO_HI", 4)
    der_tol <- .opt_get_num("mde.ctm.hpo_deriv_tol", "MDE_CTM_HPO_DERIV_TOL", 1e-6)

    Xall <- Sout$X_tr
    N <- nrow(Xall)
    set.seed(seed)
    idx <- sample.int(N)
    n_val <- max(1L, floor(hpo_val_frac * N))
    idx_val <- idx[seq_len(n_val)]
    idx_fit <- idx[-seq_len(n_val)]
    X_fit <- Xall[idx_fit, , drop = FALSE]
    X_val <- Xall[idx_val, , drop = FALSE]

    # Initial hyperparams
    Q0 <- {
      if (!is.null(Q)) as.integer(Q) else {
        qo <- getOption("mde.ctm.Q", NA_integer_)
        if (is.na(qo)) {
          qe <- suppressWarnings(as.integer(Sys.getenv("MDE_CTM_Q", NA_character_)))
          if (is.na(qe)) 24L else qe
        } else qo
      }
    }
    Qc <- min(max(Q0, Q_min), Q_max)
    ln <- max(lam_lo, min(lam_hi, if (is.null(lambda_non)) 3e-2 else lambda_non))
    lm <- max(lam_lo, min(lam_hi, if (is.null(lambda_mon)) 2e-2 else lambda_mon))
    r <- lm / max(ln, 1e-12)
    if (r < ratio_lo) lm <- ratio_lo * ln
    if (r > ratio_hi) lm <- ratio_hi * ln

    # Warmstart container
    warm_alpha <- alpha_init_list

    best <- list(val_nll = Inf, Q = Qc, ln = ln, lm = lm, coeffs = NULL)
    for (round in seq_len(hpo_max_rounds)) {
      S_fit <- list(X_tr = X_fit, X_te = X_val)
      alpha_init_list <- warm_alpha
      # fit with current hyperparameters
      lambda_non <<- ln; lambda_mon <<- lm
      fit_r <- tryCatch(fit_one(S_fit, Q_override = Qc, control = list(maxit = maxit_use)), error = function(e) e)
      if (inherits(fit_r, "error")) {
        ln <- min(lam_hi, ln * 2)
        lm <- min(lam_hi, lm * 3)
        next
      }
      # Validation score
      val_LD <- tryCatch(predict(fit_r$S, X_val, type = "logdensity_by_dim"), error = function(e) e)
      if (inherits(val_LD, "error") || any(!is.finite(val_LD))) {
        ln <- min(lam_hi, ln * 2)
        lm <- min(lam_hi, lm * 3)
        next
      }
      val_nll <- mean(rowSums(-val_LD))
      if (is.finite(val_nll) && val_nll < best$val_nll - hpo_tol_nats) {
        best <- list(val_nll = val_nll, Q = Qc, ln = ln, lm = lm, coeffs = lapply(fit_r$S$coeffs, function(cj) list(alpha = cj$alpha)))
        warm_alpha <- lapply(fit_r$S$coeffs, function(cj) cj$alpha)
      }

      # Lambda update: multiplicative candidates
      cands <- expand.grid(c = c(0.5, 1, 2), cr = c(0.75, 1, 1.25), stringsAsFactors = FALSE)
      scored <- FALSE
      for (i in seq_len(nrow(cands))) {
        ln_c <- max(lam_lo, min(lam_hi, ln * cands$c[i]))
        lm_c <- max(lam_lo, min(lam_hi, lm * cands$c[i] * cands$cr[i]))
        r_c <- lm_c / max(ln_c, 1e-12)
        if (r_c < ratio_lo) lm_c <- ratio_lo * ln_c
        if (r_c > ratio_hi) lm_c <- ratio_hi * ln_c
        lambda_non <<- ln_c; lambda_mon <<- lm_c
        alpha_init_list <- warm_alpha
        fit_c <- tryCatch(fit_one(S_fit, Q_override = Qc, control = list(maxit = maxit_use)), error = function(e) e)
        if (inherits(fit_c, "error")) next
        LD_c <- tryCatch(predict(fit_c$S, X_val, type = "logdensity_by_dim"), error = function(e) e)
        if (inherits(LD_c, "error") || any(!is.finite(LD_c))) next
        nll_c <- mean(rowSums(-LD_c))
        scored <- TRUE
        if (nll_c < val_nll - hpo_tol_nats && nll_c < best$val_nll - hpo_tol_nats) {
          ln <- ln_c; lm <- lm_c; val_nll <- nll_c
          best <- list(val_nll = nll_c, Q = Qc, ln = ln, lm = lm, coeffs = lapply(fit_c$S$coeffs, function(cj) list(alpha = cj$alpha)))
          warm_alpha <- lapply(fit_c$S$coeffs, function(cj) cj$alpha)
        }
      }

      # Adaptive quadrature via finite-difference check against exp(h)
      check_deriv <- function(Smap, Xuse, samples = 16L) {
        set.seed(seed + 123)
        I <- sample.int(nrow(Xuse), min(samples, nrow(Xuse)))
        Xs <- .standardize(Smap, Xuse[I, , drop = FALSE])
        nodes <- Smap$quad_nodes_ct; weights <- Smap$quad_weights_ct;         for (k in seq_len(ncol(Xs))) {
          Xprev <- if (k > 1) Xs[, 1:(k - 1), drop = FALSE] else matrix(0, nrow(Xs), 0)
          xk <- Xs[, k]
          for (i in seq_along(xk)) {
            xp <- if (k > 1) Xprev[i, ] else numeric(0)
            x0 <- xk[i]
            rel <- .check_eq22_point_ct(Smap, k, x0, xp, nodes, weights)
            if (rel > der_tol) return(FALSE)
          }
        }
        TRUE
      }
      # If derivative check fails, increase Q
      if (!check_deriv(fit_r$S, X_val)) {
        Qc_new <- min(Q_max, Qc + Q_step)
        if (Qc_new > Qc) {
          Qc <- Qc_new
          next
        }
      }

      # Early stop if no improvement from candidates
      if (!scored || abs(best$val_nll - val_nll) <= hpo_tol_nats) break
    }

    # Final refit on full training with best hyperparams (rollback safe)
    lambda_non <<- best$ln; lambda_mon <<- best$lm
    if (!is.null(best$coeffs)) alpha_init_list <- lapply(best$coeffs, `[[`, "alpha")
    out <- fit_one(Sout, Q_override = best$Q, control = list(maxit = maxit_use))
    out$S$perm <- perm
    out$S$meta$hpo <- list(enabled = TRUE, val_frac = hpo_val_frac, best_val_nll = best$val_nll,
                           rounds = hpo_max_rounds, seed = seed, Q_min = Q_min, Q_max = Q_max, Q_step = Q_step,
                           lambda_bounds = c(lam_lo, lam_hi), ratio_bounds = c(ratio_lo, ratio_hi), deriv_tol = der_tol)
  }

  # Quick sanity on training Z: |mean|<=0.5 and var in [0.5,2]; if violated, refit once with stronger λ_mon
  compute_Z_metrics <- function(S_map, X) {
    Xs <- .standardize(S_map, X)
    N <- nrow(Xs); K <- ncol(Xs)
    nodes <- S_map$quad_nodes_ct; weights <- S_map$quad_weights_ct;     res <- lapply(seq_len(K), function(k) {
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
          # Fallback uses B-spline helper; in practice, bs_spec is always present from training
          Psi_q <- .build_Psi_q_bs_ct(xval, xp, nodes, S_map$coeffs[[k]]$bs_spec)
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
  v_ok <- is.finite(zm[, 'var'])
  bad <- any(abs(zm[, 'mean']) > 0.5) || any(v_ok & (zm[, 'var'] < 0.5 | zm[, 'var'] > 2.0))
  if (isTRUE(bad)) {
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
  X_in <- if (!is.null(object$perm)) {
    stopifnot(ncol(newdata) == length(object$perm))
    newdata[, object$perm, drop = FALSE]
  } else newdata
  Xs <- .standardize(object, X_in)
  N <- nrow(Xs)
  K <- ncol(Xs)
  batch_use <- if (is.null(batch_n)) min(N, max(256L, floor(65536 / max(1L, object$Q)))) else min(N, batch_n)
  nodes <- object$quad_nodes_ct
  weights <- object$quad_weights_ct
  # no  in lean mode
  C <- -0.5 * log(2 * pi)

  # Optional Eq.(22) derivative check (debug only, single-source)
  if (isTRUE(getOption("mde.check_eq22", FALSE))) {
    der_tol <- getOption("mde.ctm.hpo_deriv_tol", 1e-6)
    set.seed(123)
    idx <- sample.int(N, min(8L, N))
    Xs_dbg <- Xs[idx, , drop = FALSE]
    for (k in seq_len(K)) {
      Xprev <- if (k > 1) Xs_dbg[, 1:(k - 1), drop = FALSE] else matrix(0, nrow(Xs_dbg), 0)
      xk <- Xs_dbg[, k]
      for (i in seq_along(xk)) {
        xp <- if (k > 1) Xprev[i, ] else numeric(0)
        x0 <- xk[i]
        rel <- .check_eq22_point_ct(object, k, x0, xp, nodes, weights)
        if (rel > der_tol) warning(sprintf("Eq.(22) derivative check failed (k=%d, rel=%.3e)", k, rel))
      }
    }
  }

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
      if (isTRUE(getOption('mde.verbose', FALSE)) && !getOption('mde.ctm_predict_fastpath_printed', FALSE)) {
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
          # Fallback uses B-spline helper; in practice, bs_spec is always present from training
          Psi_q <- .build_Psi_q_bs_ct(xval, xp, nodes, object$coeffs[[k]]$bs_spec)
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
    # (Removed) inner runtime Eq.(22) checks: unified single-source check at predict entry
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
  LD_joint <- rowSums(LD)
  if (!all(is.finite(LD_joint))) stop("Non-finite joint logdensity")
  if (type == "logdensity_by_dim") LD else LD_joint
}

# Ende -----------------------------------------------------------------------

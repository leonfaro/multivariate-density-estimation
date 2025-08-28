# Normalizing Flow (RealNVP-style, shallow linear couplings)
# Proof-of-concept in base R; no external deps; train/test only.

if (!exists(".log2pi")) .log2pi <- log(2 * pi)

.nf_make_masks <- function(K, L) {
  # Alternate odd/even masks per layer
  odd <- seq(1, K, by = 2)
  even <- setdiff(seq_len(K), odd)
  masks <- vector("list", L)
  for (l in seq_len(L)) {
    A <- if (l %% 2 == 1) odd else even
    B <- setdiff(seq_len(K), A)
    masks[[l]] <- list(A = A, B = B)
  }
  masks
}

.nf_param_sizes <- function(K, masks) {
  # For each layer, params for s,t: W_s (q x p), b_s(q), W_t(q x p), b_t(q)
  sizes <- lapply(masks, function(m) {
    p <- length(m$A); q <- length(m$B)
    list(p = p, q = q,
         W_s = q * p, b_s = q, W_t = q * p, b_t = q,
         total = 2 * q * p + 2 * q)
  })
  list(sizes = sizes, total = sum(vapply(sizes, `[[`, integer(1), "total")))
}

.nf_unpack_theta <- function(theta, K, masks) {
  ps <- .nf_param_sizes(K, masks)$sizes
  off <- 0L
  params <- vector("list", length(masks))
  for (l in seq_along(masks)) {
    p <- ps[[l]]$p; q <- ps[[l]]$q
    W_s <- matrix(theta[(off + 1):(off + q * p)], nrow = q, ncol = p); off <- off + q * p
    b_s <- theta[(off + 1):(off + q)]; off <- off + q
    W_t <- matrix(theta[(off + 1):(off + q * p)], nrow = q, ncol = p); off <- off + q * p
    b_t <- theta[(off + 1):(off + q)]; off <- off + q
    params[[l]] <- list(W_s = W_s, b_s = b_s, W_t = W_t, b_t = b_t)
  }
  stopifnot(off == length(theta))
  params
}

.nf_forward <- function(X, K, masks, params) {
  # Apply L coupling layers; return Z and log|det J|
  Xcur <- X
  N <- nrow(Xcur)
  Lsum <- rep(0, N)
  for (l in seq_along(masks)) {
    m <- masks[[l]]; par <- params[[l]]
    A <- m$A; B <- m$B
    XA <- Xcur[, A, drop = FALSE]; XB <- Xcur[, B, drop = FALSE]
    S <- sweep(XA %*% t(par$W_s), 2, par$b_s, "+")
    T <- sweep(XA %*% t(par$W_t), 2, par$b_t, "+")
    YB <- XB * exp(S) + T
    Xnext <- Xcur
    Xnext[, B] <- YB
    Xcur <- Xnext
    Lsum <- Lsum + rowSums(S)
  }
  list(Z = Xcur, logdet = Lsum)
}

.nf_negloglik <- function(theta, X, K, masks) {
  params <- .nf_unpack_theta(theta, K, masks)
  fw <- .nf_forward(X, K, masks, params)
  Z <- fw$Z; logdet <- fw$logdet
  base <- rowSums((-0.5) * (Z^2) - 0.5 * .log2pi)
  nll <- -mean(base + logdet)
  if (!is.finite(nll)) nll <- Inf
  nll
}

trainNFSurrogate <- function(S_or_X, layers = 2L, seed = 42) {
  set.seed(seed)
  X <- if (is.list(S_or_X)) S_or_X$X_tr else S_or_X
  stopifnot(is.matrix(X))
  K <- ncol(X)
  masks <- .nf_make_masks(K, layers)
  total <- .nf_param_sizes(K, masks)$total
  theta0 <- rnorm(total, sd = 0.05)
  fn <- function(th) .nf_negloglik(th, X, K, masks)
  opt <- optim(theta0, fn, method = "L-BFGS-B", control = list(maxit = 300))
  params <- .nf_unpack_theta(opt$par, K, masks)
  model <- list(K = K, masks = masks, params = params, seed = seed, layers = layers)
  class(model) <- "nf_surrogate"
  model
}

predict.nf_surrogate <- function(object, newdata, type = c("logdensity_by_dim", "logdensity")) {
  type <- match.arg(type)
  X <- as.matrix(newdata)
  fw <- .nf_forward(X, object$K, object$masks, object$params)
  Z <- fw$Z; logdet <- fw$logdet
  base_by_dim <- (-0.5) * (Z^2) - 0.5 * .log2pi  # N x K
  # Distribute logdet to scaled dims per layer (accumulate contributions per dim)
  N <- nrow(X); K <- ncol(X)
  logdet_by_dim <- matrix(0, N, K)
  Xcur <- X
  for (l in seq_along(object$masks)) {
    m <- object$masks[[l]]; par <- object$params[[l]]
    A <- m$A; B <- m$B
    XA <- Xcur[, A, drop = FALSE]
    S <- sweep(XA %*% t(par$W_s), 2, par$b_s, "+")
    logdet_by_dim[, B] <- logdet_by_dim[, B, drop = FALSE] + S
    # update Xcur like in forward to keep consistency for next layer
    XB <- Xcur[, B, drop = FALSE]
    T <- sweep(XA %*% t(par$W_t), 2, par$b_t, "+")
    YB <- XB * exp(S) + T
    Xcur[, B] <- YB
  }
  LD <- base_by_dim + logdet_by_dim
  stopifnot(all(is.finite(LD)))
  if (type == "logdensity_by_dim") return(LD)
  rowSums(LD)
}

# Inverse pass and sampling --------------------------------------------------

.nf_inverse <- function(Z, K, masks, params) {
  Ycur <- Z
  for (l in rev(seq_along(masks))) {
    m <- masks[[l]]; par <- params[[l]]
    A <- m$A; B <- m$B
    YA <- Ycur[, A, drop = FALSE]
    YB <- Ycur[, B, drop = FALSE]
    S <- sweep(YA %*% t(par$W_s), 2, par$b_s, "+")
    T <- sweep(YA %*% t(par$W_t), 2, par$b_t, "+")
    XB <- (YB - T) * exp(-S)
    Xprev <- Ycur
    Xprev[, B] <- XB
    Ycur <- Xprev
  }
  Ycur
}

sample.nf_surrogate <- function(object, n) {
  Z <- matrix(rnorm(n * object$K), ncol = object$K)
  .nf_inverse(Z, object$K, object$masks, object$params)
}

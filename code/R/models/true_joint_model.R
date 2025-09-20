# TRUE joint model: oracle evaluation of triangular joint density

# sanitize parameter list for distribution
.sanitize_args <- function(args, distr) {
  if (distr == "gamma" && all(c("shape1", "shape2") %in% names(args))) {
    args$shape <- args$shape1
    args$scale <- args$shape2
    args$shape1 <- NULL
    args$shape2 <- NULL
  }
  pos_params <- switch(distr,
    norm  = "sd",
    exp   = "rate",
    beta  = c("shape1", "shape2"),
    gamma = c("shape", "scale"),
    character(0)
  )
  defaults <- list(
    norm = list(sd = 1),
    exp = list(rate = 1),
    beta = list(shape1 = 1, shape2 = 1),
    gamma = list(shape = 1, scale = 1)
  )
  for (nm in pos_params) {
    val <- args[[nm]]
    if (is.null(val) || !is.finite(val) || val <= 0) {
      val <- defaults[[distr]][[nm]]
    }
    args[[nm]] <- max(val, 1e-6)
  }
  args
}

# compute conditional log-density for single observation row
.log_density_conditional_row <- function(x_row, config) {
  K <- length(config)
  out <- numeric(K)
  prev_names <- paste0("X", seq_len(K))
  for (k in seq_len(K)) {
    prev <- if (k == 1) data.frame() else {
      df <- as.data.frame(as.list(x_row[seq_len(k - 1)]))
      names(df) <- prev_names[seq_len(k - 1)]
      df
    }
    distr_k <- config[[k]]$distr
    if (is.null(config[[k]]$parm)) {
      args <- list()
      if (distr_k == "beta") {
        args$shape1 <- 1
        args$shape2 <- 1
      } else if (distr_k == "gamma") {
        args$shape <- 1
        args$scale <- 1
      }
    } else {
      args <- config[[k]]$parm(prev)
    }
    args <- .sanitize_args(args, distr_k)
    xk <- x_row[k]
    if (distr_k %in% c("exp", "gamma")) {
      xk <- max(xk, 1e-6)
    } else if (distr_k == "beta") {
      xk <- min(max(xk, 1e-6), 1 - 1e-6)
    }
    out[k] <- switch(distr_k,
      norm  = dnorm(xk, mean = args$mean %||% 0, sd = args$sd %||% 1, log = TRUE),
      exp   = dexp(xk, rate = args$rate %||% 1, log = TRUE),
      beta  = dbeta(xk, shape1 = args$shape1, shape2 = args$shape2, log = TRUE),
      gamma = dgamma(xk, shape = args$shape, scale = args$scale, log = TRUE),
      stop("Unsupported distribution")
    )
  }
  out
}

#' Row-wise conditional log-densities
#' @param config configuration list
#' @param X matrix of observations
#' @param cores parallel cores
#' @return N x K matrix of log-densities
#' @export
true_joint_logdensity_by_dim <- function(config, X, cores = NC) {
  stopifnot(is.matrix(X))
  res <- parallel::mclapply(seq_len(nrow(X)), function(i) {
    .log_density_conditional_row(X[i, ], config)
  }, mc.cores = cores)
  ll <- do.call(rbind, res)
  if (!all(is.finite(ll))) stop("log-density not finite")
  ll
}

#' Dimension-wise negative log-likelihood under true joint
#' @export
logL_TRUE_JOINT_dim <- function(config, X, cores = NC) {
  ll <- true_joint_logdensity_by_dim(config, X, cores)
  -colMeans(ll)
}

#' Overall negative log-likelihood under true joint
#' @export
logL_TRUE_JOINT <- function(config, X, cores = NC) {
  ll <- true_joint_logdensity_by_dim(config, X, cores)
  -mean(rowSums(ll))
}

#' Evaluate true joint model on test data
#' @param S list with element X_te
#' @param config configuration list
#' @return list with class 'true_joint'
#' @export
fit_TRUE_JOINT <- function(S, config, cores = NC) {
  stopifnot(is.list(S))
  te_dim <- logL_TRUE_JOINT_dim(config, S$X_te, cores)
  res <- list(
    config = config,
    logL_te_dim = te_dim,
    logL_te = sum(te_dim)
  )
  class(res) <- "true_joint"
  res
}

`%||%` <- function(a, b) if (is.null(a)) b else a


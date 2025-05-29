#' Generate conditional samples
#'
#' This function implements Algorithm 1 as described in `roadmap.md`.
#' It sequentially draws samples from the distributions specified in
#' `config`, potentially conditioning on previously generated columns.
#'
#' @param G list with elements `N`, `config` and `seed`
#' @return numeric matrix of dimension `N` x `length(config)`
#' @examples
#' G <- setup_global()
#' X <- gen_samples(G)
#' dim(X)
#'

# helper to draw one observation from distribution `distr`
.draw_from <- function(distr, params) {
  fun <- get(paste0("r", distr), mode = "function")
  if (distr == "gamma" && all(c("shape1", "shape2") %in% names(params))) {
    params <- list(shape = params$shape1, scale = params$shape2)
  }
  params <- lapply(params, function(p) ifelse(p <= 0, 1e-3, p))
  do.call(fun, c(list(n = 1), params))
}

#' @rdname gen_samples
#' @export
gen_samples <- function(G) {
  N <- G$N
  config <- G$config
  seed <- G$seed
  K <- length(config)

  set.seed(seed)
  X <- matrix(NA_real_, nrow = N, ncol = K)

  for (i in seq_len(N)) {
    for (k in seq_len(K)) {
      c_k <- config[[k]]
      if (is.null(c_k$parm)) {
        params_k <- list()
      } else {
        if (k == 1) {
          prev <- data.frame()
        } else {
          prev <- as.data.frame(as.list(X[i, seq_len(k - 1)]))
          names(prev) <- paste0("X", seq_len(k - 1))
        }
        params_k <- c_k$parm(prev)
      }
      X[i, k] <- .draw_from(c_k$distr, params_k)
    }
  }

  colnames(X) <- paste0("X", seq_len(K))
  X
}

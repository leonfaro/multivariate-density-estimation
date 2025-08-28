#' Generate conditional samples
#'
#' This function implements conditional sampling via triangular transport mapping as described in `Theory.md`.
#' It sequentially draws samples from the distributions specified in
#' `config`, potentially conditioning on previously generated columns.
#'

# 0.1  Bibliotheks-freie Pseudocode-Hilfsfunktionen -------------------

Generate_iid_from_config <- function(N, cfg, return_params = FALSE) {
  stopifnot(is.numeric(N), N > 0, is.list(cfg))
  K <- length(cfg)
  X <- matrix(NA_real_, nrow = N, ncol = K)
  if (return_params) {
    param_hist <- vector("list", K)
    for (kk in seq_len(K))
      param_hist[[kk]] <- vector("list", N)
  }
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
      if (return_params) param_hist[[k]][[i]] <- args
      X[i, k] <- do.call(fun, c(list(n = 1L), args))
    }
  }
  colnames(X) <- paste0("X", seq_len(K))
  if (return_params) {
    param_df <- lapply(param_hist, function(lst) {
      vals <- lapply(lst, function(x) if (length(x) == 0) NULL else as.data.frame(x))
      vals <- Filter(Negate(is.null), vals)
      if (length(vals) == 0) return(NULL)
      df <- do.call(rbind, vals)
      rownames(df) <- NULL
      df
    })
    list(X = X, params = param_df)
  } else {
    X
  }
}


# Kompatibilitäts-Funktion für bestehende Tests
gen_samples <- function(G, return_params = FALSE) {
  Generate_iid_from_config(G$n, G$config, return_params = return_params)
}

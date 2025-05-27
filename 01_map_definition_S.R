# Map definition functions
dist_fun <- function(pref, name) get(paste0(pref, name))
get_pars <- function(k, x_prev, cfg) {
  ck <- cfg[[k]]
  if (is.null(ck$parm)) return(list())
  if (is.null(dim(x_prev))) {
    names(x_prev) <- paste0("X", seq_along(x_prev))
    x_df <- as.data.frame(as.list(x_prev))
  } else {
    colnames(x_prev) <- paste0("X", seq_len(ncol(x_prev)))
    x_df <- as.data.frame(x_prev)
  }
  pars <- ck$parm(x_df)
  pars <- apply_links(pars, ck$distr)
  safe_pars(pars, ck$distr)
}

pdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  pars <- get_pars(k, x_prev, cfg)
  xk   <- safe_support(xk, cfg[[k]]$distr, pars)
  dens <- do.call(dist_fun("d", cfg[[k]]$distr),
                  c(list(x = xk), pars, list(log = log)))
  dens
}

cdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  dname <- cfg[[k]]$distr
  pars  <- get_pars(k, x_prev, cfg)
  xk    <- safe_support(xk, dname, pars)
  support <- p_supports_logp[dname]
  if (is.na(support))
    stop("log.p capability not specified for distribution ", dname)
  if (log && support) {
    args <- c(list(q = xk),
              pars,
              list(log.p = TRUE))
    cdfv <- do.call(dist_fun("p", dname), args)
  } else {
    args <- c(list(q = xk),
              pars,
              list(log.p = FALSE))
    cdfv <- do.call(dist_fun("p", dname), args)
    if (log) log(cdfv) else cdfv
  }
}

qtf_k <- function(k, u, x_prev, cfg, log.p = FALSE) {
  dname <- cfg[[k]]$distr
  pars <- get_pars(k, x_prev, cfg)
  support <- q_supports_logp[dname]
  if (is.na(support))
    stop("log.p capability not specified for distribution ", dname)
  if (log.p && !support) {
    u <- exp(u)
    log.p <- FALSE
  }
  res <- do.call(dist_fun("q", dname), c(list(p = u, log.p = log.p), pars))
  safe_support(res, dname, pars)
}

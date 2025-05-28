# Input: dist_registry, link_fns, get_pars(), config
# Output: pdf_k(), cdf_k(), qtf_k()
# Map definition functions

pdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  pars <- get_pars(k, x_prev, cfg)
  dens <- do.call(dist_fun("d", cfg[[k]]$distr),
                  c(list(x = xk), pars, list(log = log)))
  dens
}

cdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  dname <- cfg[[k]]$distr
  pars  <- get_pars(k, x_prev, cfg)
  do.call(dist_fun("p", dname),
          c(list(q = xk, log.p = log), pars))
}

qtf_k <- function(k, u, x_prev, cfg, log.p = FALSE) {
  dname <- cfg[[k]]$distr
  pars <- get_pars(k, x_prev, cfg)
  do.call(dist_fun("q", dname), c(list(p = u, log.p = log.p), pars))
}

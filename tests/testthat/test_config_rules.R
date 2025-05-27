# Übersicht – Differenzierbarkeit & Träger für 18 stetige Verteilungen
# Legende:
# glatt = CInf, cusp = Spitze/Knick (nicht C1)
#
# 1 Normal:         (-Inf,Inf); mu in R, sigma>0; glatt; glatt
# 2 Lognormal:      (0,Inf);   mu in R, sigma>0; cusp bei 0+; glatt innen
# 3 Student-t:      (-Inf,Inf); mu in R, sigma>0, nu>0; glatt; glatt
# 4 Skew-t Hansen:  (-Inf,Inf); mu in R, sigma>0, nu>2, lambda in(-1,1); glatt
# 5 GED:            (-Inf,Inf); mu in R, beta>0, p>0; p<=1 cusp bei mu
# 6 NIG:            (-Inf,Inf); alpha>0, |beta|<alpha, delta>0, mu in R; glatt
# 7 Variance-Gamma: (-Inf,Inf); kappa>0, theta in R, sigma>0, mu in R; glatt
# 8 alpha-Stable:   (-Inf,Inf); alpha in(0,2], beta in[-1,1], gamma>0, mu in R;
#                    alpha>1 glatt, sonst cusp bei mu
# 9 Exponential:    (0,Inf); rate>0; glatt auch am Rand
#10 Weibull:        (0,Inf); k>0, lambda>0; k<1 cusp bei 0+
#11 Laplace:        (-Inf,Inf); mu in R, b>0; cusp bei mu
#12 Gen. Pareto:    [mu, mu - sigma/xi) falls xi<0 sonst [mu, Inf); sigma>0;
#                    xi in R; xi<1 cusp bei mu
#13 Burr XII:       (0,Inf); c>0, k>0, lambda>0; c<1 cusp bei 0+
#14 Hyperbolic:     (-Inf,Inf); alpha>0, |beta|<alpha, delta>0, mu in R; glatt
#15 Inverse-Gauss.: (0,Inf); mu>0, lambda>0; cusp bei 0+
#16 Non-central chi^2: (0,Inf); nu>0, lambda>=0; nu<2 cusp bei 0+
#17 Gamma:          (0,Inf); alpha>0, theta>0; alpha<1 cusp bei 0+
#18 Beta:           (0,1); alpha>0, beta>0; glatt wenn alpha>1 & beta>1

library(testthat)

extract_config <- function(file) {
  exprs <- parse(file)
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name("config"))) {
      return(eval(e[[3]], envir = new.env()))
    }
  }
  stop("config not found in ", file)
}

check_cfg <- function(cfg, root) {
  env <- new.env()
  env$config <- cfg
  sys.source(file.path(root, "00_setup.R"), env, chdir = TRUE)
  allowed_dists <- c(
    "norm", "lnorm", "t", "skewt", "ged", "nig", "vg",
    "astable", "exp", "weibull", "laplace", "gpd", "burrXII",
    "hyperbolic", "invgauss", "ncchisq", "gamma", "beta"
  )
  for (k in seq_along(cfg)) {
    ck <- cfg[[k]]
    expect_true(is.list(ck))
    expect_true(is.character(ck$distr))
    expect_true(ck$distr %in% names(env$q_supports_logp))
    expect_true(ck$distr %in% allowed_dists)

    if (!is.null(ck$parm)) {
      fn_txt <- paste(deparse(ck$parm), collapse = " ")
      # Check no future dimensions referenced
      max_dim <- max(length(cfg), 6)
      for (d in seq(from = k, to = max_dim)) {
        expect_false(grepl(paste0("X", d), fn_txt, fixed = TRUE))
      }
      # Positive parameter enforcement via monotone transform
      need_pos <- ck$distr %in% c(
        "norm", "exp", "gamma", "weibull", "lnorm", "pois", "beta",
        "logis", "laplace", "t", "skewt", "ged", "nig", "vg", "astable",
        "gpd", "burrXII", "hyperbolic", "invgauss", "ncchisq"
      )
      if (need_pos) {
        expect_true(grepl("softplus", fn_txt) ||
                    grepl("exp(", fn_txt, fixed = TRUE) ||
                    grepl("plogis", fn_txt))
      }
      # Call function on zero data frame
      if (k > 1) {
        d <- setNames(as.data.frame(matrix(0, nrow = 1, ncol = k - 1)),
                       paste0("X", seq_len(k - 1)))
      } else {
        d <- data.frame()[, FALSE]
      }
      environment(ck$parm) <- env
      pars <- ck$parm(d)
      expect_silent(env$safe_pars(pars, ck$distr))
    }
  }
}

test_that("configs satisfy monotonicity rules", {
  root <- file.path("..", "..")
  check_cfg(extract_config(file.path(root, "run3.R")), root)
  check_cfg(extract_config(file.path(root, "run5.R")), root)
  check_cfg(extract_config(file.path(root, "basic_tests.R")), root)
})

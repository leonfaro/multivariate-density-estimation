library(testthat)

# Übersicht – Differenzierbarkeit & Träger für die 17 relevanten stetigen
# Verteilungen
#
# Legende
# * glatt   = Klasse C∞ (beliebig oft stetig differenzierbar)
# * C¹ / C² = nur bis zur genannten Ordnung stetig differenzierbar
# * cusp    = Spitze oder Knick, Ableitungen nicht endlich
# "Rand"    = linkes bzw. rechtes Intervallende des Supports
#
# | #  | Verteilung       | a) Support                                                 | b) Parameter­bereich                           | c) Glattheit auf ganzem Support                                | d) Glattheit im Inneren             |
# |----|------------------|-----------------------------------------------------------|------------------------------------------------|---------------------------------------------------------------|------------------------------------|
# | 1  | Normal           | (-∞,∞)                                                   | μ∈ℝ, σ>0                                       | glatt                                                         | glatt                              |
# | 2  | Lognormal        | (0,∞)                                                   | μ∈ℝ, σ>0                                       | cusp am linken Rand x→0⁺ (nicht C¹)                          | glatt                              |
# | 3  | Student-t        | (-∞,∞)                                                   | μ∈ℝ, σ>0, ν>0                                  | glatt                                                         | glatt                              |
# | 4  | Skew-t (Hansen)  | (-∞,∞)                                                   | μ∈ℝ, σ>0, ν>2, λ∈(-1,1)                        | glatt                                                         | glatt                              |
# | 5  | GED              | (-∞,∞)                                                   | μ∈ℝ, β>0, p>0                                  | p>1: glatt; p=1: cusp bei x=μ; p<1: cusp stärker              | analog, glatt nur falls p>1        |
# | 6  | NIG              | (-∞,∞)                                                   | α>0, β<α, δ>0, μ∈ℝ                             | glatt                                                         | glatt                              |
# | 7  | Variance-Gamma   | (-∞,∞)                                                   | κ>0, θ∈ℝ, σ>0, μ∈ℝ                             | glatt                                                         | glatt                              |
# | 8  | α-Stable         | (-∞,∞)                                                   | α∈(0,2], β∈[-1,1], γ>0, μ∈ℝ                    | α>1: C¹ (praktisch glatt); α≤1: cusp bei x=μ                  | α>1: glatt; α≤1: cusp              |
# | 10a| Exponential      | (0,∞)                                                    | λ>0                                           | glatt (auch am Rand)                                         | glatt                              |
# | 10b| Weibull          | (0,∞)                                                    | k>0, λ>0                                       | k≥1: glatt; k<1: cusp bei 0⁺                                   | k≥1: glatt; k<1: glatt ab x>0      |
# | 11 | Laplace          | (-∞,∞)                                                   | μ∈ℝ, b>0                                       | cusp bei x=μ (nicht C¹)                                      | cusp bleibt                       |
# | 12 | Gen. Pareto      | [μ, μ-σ/ξ) falls ξ<0 sonst [μ,∞)                          | σ>0, ξ∈ℝ, μ∈ℝ                                  | ξ≥1: C¹/C² je nach ξ; ξ<1: cusp bei x=μ                       | glatt für x>μ bzw. <Obergrenze     |
# | 13 | Burr XII         | (0,∞)                                                    | c>0, k>0, λ>0                                  | c≥1: glatt; c<1: cusp bei 0⁺                                  | glatt ab x>0                       |
# | 14 | Hyperbolic       | (-∞,∞)                                                   | α>0, β<α, δ>0, μ∈ℝ                             | glatt                                                         | glatt                              |
# | 15 | Inverse-Gaussian | (0,∞)                                                    | μ>0, λ>0                                       | cusp bei 0⁺ (x^{-3/2}-Singularität)                          | glatt für x>0                      |
# | 16 | Non-central χ²   | (0,∞)                                                    | ν>0, λ≥0                                       | ν≥2: glatt; ν<2: cusp bei 0⁺                                   | glatt für x>0                      |
# | 17 | Gamma            | (0,∞)                                                    | α>0, θ>0                                       | α≥1: C¹ (für α≥2 glatt); α<1: cusp bei 0⁺                     | glatt ab x>0                       |
# | 18 | Beta             | (0,1)                                                    | α>0, β>0                                       | glatt falls α>1 & β>1; sonst cusp an 0 bzw. 1                 | glatt auf (0,1)                    |
#
# Hinweise zu Randproblemen:
# "cusp" signalisiert eine nicht endliche Ableitung am Rand oder an einem
# inneren Knickpunkt. Glatte Optimierungsverfahren benötigen daher passende
# Reparametrisierungen oder Ausschlüsse solcher Parameterbereiche.

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
    "norm", "lnorm", "t", "skewt",
    "exp", "weibull", "laplace", "gpd", "burrxii",
    "invgauss", "ncchisq", "gamma", "beta", "logis"
  )
  pos_pars <- list(
    norm    = "sigma",
    logis   = "scale",
    t       = c("df", "sigma"),
    cauchy  = "scale",
    gumbel  = "scale",
    lnorm   = "sdlog",
    gamma   = c("shape", "rate"),
    weibull = c("shape", "scale"),
    beta    = c("shape1", "shape2"),
    exp     = "rate"
  )
  support_map <- c(
    norm = "real", logis = "real", t = "real", cauchy = "real", gumbel = "real",
    lnorm = "positive", gamma = "positive", weibull = "positive",
    beta = "unit", exp = "positive0"
  )
  sup_samples <- list(
    real = c(-1, 0, 1),
    positive = c(0.1, 1, 2),
    positive0 = c(0, 1, 2),
    unit = c(0.1, 0.5, 0.9)
  )


  for (k in seq_along(cfg)) {
    ck <- cfg[[k]]
    expect_true(is.list(ck))
    expect_true(is.character(ck$distr))

    expect_true(ck$distr %in% allowed_dists)
    if (!is.null(ck$parm)) {
      fn_txt <- paste(deparse(ck$parm), collapse = " ")
      # Check no future dimensions referenced
      max_dim <- max(length(cfg), 6)
      for (d in seq(from = k, to = max_dim)) {
        expect_false(grepl(paste0("X", d), fn_txt, fixed = TRUE))
      }
      need_pos <- ck$distr %in% names(pos_pars)
      # Call function on zero data frame
      if (k > 1) {
        d <- setNames(as.data.frame(matrix(0, nrow = 1, ncol = k - 1)),
                       paste0("X", seq_len(k - 1)))
      } else {
        d <- data.frame()[, FALSE]
      }
      environment(ck$parm) <- env
      x_prev <- if (k > 1) as.numeric(d) else numeric(0)
      expect_silent(pars <- env$get_pars(k, x_prev, cfg))

      if (need_pos) {
        for (p in pos_pars[[ck$distr]]) {
          if (!is.null(pars[[p]]))
            expect_true(all(pars[[p]] > 0))
        }
      }
      if (ck$distr %in% names(env$dist_registry)) {
        sp <- support_map[[ck$distr]]
        xs <- sup_samples[[sp]]
        logpdf_fun <- env$dist_registry[[ck$distr]]$logpdf
        vals <- do.call(logpdf_fun, c(list(x = xs), pars))
        expect_true(all(is.finite(vals)))
      }


    }
  }
}

test_that("configs satisfy monotonicity rules", {
  root <- file.path("..", "..")
  check_cfg(extract_config(file.path(root, "run3.R")), root)
  run5 <- file.path(root, "run5.R")
  if (file.exists(run5)) {
    check_cfg(extract_config(run5), root)
  } else {
    skip("run5.R not available")
  }
  check_cfg(extract_config(file.path(root, "basic_tests.R")), root)
})

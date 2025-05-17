## ─────────────────────────  Conditional-chain config25  ─────────────────────────
##  • X1  ~ f1(·)  (unconditional anchor)
##  • Xk ~ fk(· | X{k-1})  for k = 2,…,25    ← strict one-step linear dependence
##  • Parameters clipped in (lo, hi) to avoid zeros / Infs
##  • All distributions keep the original ordering

library(extraDistr)   # dpower, dlomax, etc.
library(sn)           # dsn / qsn   (for k = 16)


config25 <- list(
  ## 1  Normal  – start node
  list(distr = "norm",
       parm  = NULL),                                   # mean = 0, sd = 1

  ## 2  Student-t   df = 3 + 0.5·X1
  list(distr = "t",
       parm  = \(d) list(df = 3 + 0.5*d$X1)),

  ## 3  Laplace     μ = 0.3·X2, σ = 1 + 0.1·X2
  list(distr = "laplace",
       parm  = \(d) list(mu    = 0.3*d$X2,
                         sigma = 1 + 0.1*d$X2)),

  ## 4  Logistic    loc = 0.2·X3, scale = 1 + 0.05·X3
  list(distr = "logis",
       parm  = \(d) list(location = 0.2*d$X3,
                         scale    = 1 + 0.05*d$X3)),

  ## 5  Cauchy      loc = 0.1·X4, scale = 0.5 + 0.05·X4
  list(distr = "cauchy",
       parm  = \(d) list(location = 0.1*d$X4,
                         scale    = 0.5 + 0.05*d$X4)),

  ## 6  Exponential λ = 1 + 0.1·X5
  list(distr = "exp",
       parm  = \(d) list(rate = 1 + 0.1*d$X5)),

  ## 7  Gamma       shape = 2 + 0.2·X6,  rate = 1
  list(distr = "gamma",
       parm  = \(d) list(shape = 2 + 0.2*d$X6,
                         rate  = 1)),

  ## 8  Weibull     shape = 2 + 0.2·X7,  scale = 1 + 0.1·X7
  list(distr = "weibull",
       parm  = \(d) list(shape = 2 + 0.2*d$X7,
                         scale = 1 + 0.1*d$X7)),

  ## 9  Log-Normal  μlog = 0.3·X8,  σlog = 0.5 + 0.05·X8
  list(distr = "lnorm",
       parm  = \(d) list(meanlog = 0.3*d$X8,
                         sdlog   = 0.5 + 0.05*d$X8)),

  ## 10 Chi-square  df = 4 + 0.2·X9
  list(distr = "chisq",
       parm  = \(d) list(df = 4 + 0.2*d$X9)),

  ## 11 F           d1 = 5,   d2 = 6 + 0.3·X10
  list(distr = "f",
       parm  = \(d) list(df1 = 5,
                         df2 = 6 + 0.3*d$X10)),

  ## 12 Beta        α = 2 + 0.1·X11,  β = 2 + 0.1·X11
  list(distr = "beta",
       parm  = \(d) { s <- 2 + 0.1*d$X11
                     list(shape1 = s, shape2 = s) }),

  ## 13 Beta (skew) α = 0.5 + 0.05·X12,  β = 5 + 0.2·X12
  list(distr = "beta",
       parm  = \(d) list(shape1 = 0.5 + 0.05*d$X12,
                         shape2 = 5   + 0.2*d$X12)),

  ## 14 Gumbel      μ = 0.5·X13,  σ = 1 + 0.05·X13
  list(distr = "gumbel",
       parm  = \(d) list(mu = 0.5*d$X13,
                         sigma = 1 + 0.05*d$X13)),

  ## 15 Gen. Pareto σ = 1 + 0.05·X14, ξ = 0.1
  list(distr = "gpd",
       parm  = \(d) list(mu = 0,
                         sigma = 1 + 0.05*d$X14,
                         xi = 0.1)),

  ## 16 Skew-Normal ξ = 0.2·X15,  ω = 1,  α = 2
  list(distr = "sn",
       parm  = \(d) list(xi = 0.2*d$X15, omega = 1, alpha = 2)),

  ## 17 Power       α = 1.5 + 0.1·X16,  β = 1
  list(distr = "power",
       parm  = \(d) list(alpha = 1.5 + 0.1*d$X16,
                         beta  = 1)),

  ## 18 Inverse-Gamma α = 3 + 0.1·X17,  β = 1
  list(distr = "invgamma",
       parm  = \(d) list(alpha = 3 + 0.1*d$X17,
                         beta  = 1)),

  ## 19 Lomax       λ = 1 + 0.05·X18,  κ = 2
  list(distr = "lomax",
       parm  = \(d) list(lambda = 1 + 0.05*d$X18,
                         kappa  = 2)),

  ## 20 Rayleigh    σ = 1 + 0.05·X19
  list(distr = "rayleigh",
       parm  = \(d) list(sigma = 1 + 0.05*d$X19)),

  ## 21 Trunc-Normal μ = 0.1·X20, σ = 1
  list(distr = "tnorm",
       parm  = \(d) list(mean = 0.1*d$X20, sd = 1, a = 0, b = Inf)),

  ## 22 Frechet     λ = 1 + 0.05·X21
  list(distr = "frechet",
       parm  = \(d) list(lambda = 1 + 0.05*d$X21,
                         mu = 0, sigma = 1)),

  ## 23 Hyperb.-Secant σ = 1 + 0.05·X22
  list(distr = "hnorm",
       parm  = \(d) list(sigma = 1 + 0.05*d$X22)),

  ## 24 Normal      μ = 0.2·X23, σ = 1
  list(distr = "norm",
       parm  = \(d) list(mean = 0.2*d$X23, sd = 1)),

  ## 25 Normal      μ = 0.2·X24, σ = 1
  list(distr = "norm",
       parm  = \(d) list(mean = 0.2*d$X24, sd = 1))
)


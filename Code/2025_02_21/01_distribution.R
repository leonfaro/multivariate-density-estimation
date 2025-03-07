# 01_distribution.R

# (Neu angelegtes Skript für die Verteilungs-Definitionen)

# Modelle, jeweils $rDist und $dDist pro Dimension
#   1) Y1 ~ Chi^2(3)
#   2) Y2 ~ Gamma(shape, scale)
#   3) Y3 ~ Beta(alpha, beta)

models <- list(
  # Y1: Chi^2(3), marginal
  list(
    rDist = function(n) rchisq(n, df=3),
    dDist = function(x) dchisq(x, df=3)
  ),
  # Y2: Gamma
  list(
    rDist = function(n, shape, scale) rgamma(n, shape=shape, scale=scale),
    dDist = function(x, shape, scale) dgamma(x, shape=shape, scale=scale)
  ),
  # Y3: Beta
  list(
    rDist = function(n, alpha, beta) rbeta(n, alpha, beta),
    dDist = function(x, alpha, beta) dbeta(x, alpha, beta)
  )
)

# Bedingte Parameter je Dimension:
#   1) (NULL): Y1 ~ Chi^2(3) (marginal)
#   2) (Y1) => Y2 ~ Gamma(shape, scale)
#   3) (Y1, Y2) => Y3 ~ Beta(alpha, beta)
# Formeln so gewählt, dass Y2>0 und Y3 in [0,1].
cond <- list(
  NULL,
  function(Y1) {
    # shape = 2 + 0.1*Y1
    # scale = 0.5 + 0.01*Y1
    shp <- 2 + 0.1 * Y1
    scl <- 0.5 + 0.01 * Y1
    if(shp <= 0) shp <- 0.1
    if(scl <= 0) scl <- 0.1
    list(shape = shp, scale = scl)
  },
  function(Y1, Y2) {
    # alpha = 1 + 0.05*Y1 + 0.02*Y2
    # beta  = 1 + 0.02*Y1 + 0.02*Y2
    a_ <- 1 + 0.05*Y1 + 0.02*Y2
    b_ <- 1 + 0.02*Y1 + 0.02*Y2
    if(a_ <= 0) a_ <- 0.1
    if(b_ <= 0) b_ <- 0.1
    list(alpha = a_, beta = b_)
  }
)



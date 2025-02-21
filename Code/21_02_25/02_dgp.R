# 02_dgp.R

# Modelle, jeweils $rDist und $dDist pro Dimension
#  Hier:
#   * Y1 ~ Chi^2(3)
#   * Y2 ~ Gamma( shape, scale )
#   * Y3 ~ Beta( alpha, beta )
#
#  -> Y2 > 0 ist garantiert (Gamma) => Log(Y2) bleibt definiert
#  -> Y3 in [0,1] => Beta-Distribution

models <- list(
  # Y1: Chi^2(3), marginal
  list(
    rDist = function(n) rchisq(n, df=3),
    dDist = function(x) dchisq(x, df=3)
  ),
  # Y2: Gamma
  list(
    # shape>0, scale>0
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
# Diese Formeln wählen wir so, dass Y2 immer > 0 ist und Y3 in [0,1].
cond <- list(
  NULL,
  # Y2 param: shape = 2 + 0.1*Y1, scale = 0.5 + 0.01*Y1
  function(Y1) {
    shp <- 2 + 0.1 * Y1
    scl <- 0.5 + 0.01 * Y1
    # Minimale Sicherheitschecks
    if(shp <= 0) shp <- 0.1
    if(scl <= 0) scl <- 0.1
    list(shape = shp, scale = scl)
  },
  # Y3 param: alpha = 1 + 0.05*Y1 + 0.02*Y2, beta = 1 + 0.02*Y1 + 0.02*Y2
  function(Y1, Y2) {
    a_ <- 1 + 0.05*Y1 + 0.02*Y2
    b_ <- 1 + 0.02*Y1 + 0.02*Y2
    if(a_ <= 0) a_ <- 0.1
    if(b_ <= 0) b_ <- 0.1
    list(alpha = a_, beta = b_)
  }
)

# Generisches DGP mit Debugging-Ausgaben
my_dgp_d <- function(N=50, d=length(models),
                     models, cond,
                     chunk_size=NULL,
                     parallel=FALSE,
                     debug=FALSE) {
  
  if(debug) cat("my_dgp_d: Start with N=",N," d=",d,"\n")
  Y <- matrix(NA, nrow=N, ncol=d,
              dimnames=list(NULL, paste0("Y",1:d)))
  
  if(parallel) {
    ncores <- 2
    if(debug) cat(" -> Setting up parallel backend with",ncores,"cores\n")
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }
  
  for(j in seq_len(d)) {
    if(debug) cat("my_dgp_d: Generiere Komponente j=", j, "\n")
    
    if(j == 1) {
      # Marginale Ziehung
      if(debug) cat("  -> Marginale Ziehung aus models[[1]] (Chi^2)\n")
      Y[,1] <- models[[1]]$rDist(N)
    } else {
      # bedingte Parameter
      if(is.null(chunk_size)) {
        if(!parallel) {
          for(i in seq_len(N)) {
            paramList <- do.call(cond[[j]], as.list(Y[i,1:(j-1)]))
            Y[i,j] <- do.call(models[[j]]$rDist,
                              c(list(n=1), paramList))
          }
        } else {
          if(debug) cat("  -> Parallel loop for j=",j,"\n")
          Ycol_j <- foreach(i=1:N, .combine='c') %dopar% {
            paramList <- do.call(cond[[j]], as.list(Y[i,1:(j-1)]))
            do.call(models[[j]]$rDist, c(list(n=1), paramList))
          }
          Y[,j] <- Ycol_j
        }
      } else {
        nChunks <- ceiling(N / chunk_size)
        idx_seq <- seq_len(N)
        startIdx <- 1
        for(chunkId in seq_len(nChunks)) {
          endIdx   <- min(startIdx + chunk_size - 1, N)
          idxBlock <- idx_seq[startIdx:endIdx]
          if(debug) cat("  -> Chunk", chunkId, "Index:", startIdx, "-", endIdx, "\n")
          for(i in idxBlock) {
            paramList <- do.call(cond[[j]], as.list(Y[i,1:(j-1)]))
            Y[i,j] <- do.call(models[[j]]$rDist,
                              c(list(n=1), paramList))
          }
          startIdx <- endIdx + 1
        }
      }
    }
  }
  
  # Debug: check Y2 <= 0
  if(d>=2 && any(Y[,2] <= 0, na.rm=TRUE)) {
    cat("WARNING: In my_dgp_d, Y2 has non-positive values!\n")
    idx_bad <- which(Y[,2] <= 0)
    print(idx_bad)
    print(Y[idx_bad, 2])
  }
  
  if(debug) cat("my_dgp_d: Done, returning data.frame(Y)\n")
  data.frame(Y)
}

calc_true_dens <- function(df, models, cond) {
  n <- nrow(df)
  d <- ncol(df)
  out <- numeric(n)
  for(i in seq_len(n)) {
    prob_val <- 1
    for(j in seq_len(d)) {
      if(j == 1) {
        val <- pmax(1e-15, models[[1]]$dDist(df[i,j]))
        prob_val <- prob_val * val
      } else {
        paramList <- do.call(cond[[j]], as.list(df[i,1:(j-1)]))
        val <- pmax(1e-15, do.call(models[[j]]$dDist,
                                   c(list(df[i,j]), paramList)))
        prob_val <- prob_val * val
      }
    }
    out[i] <- prob_val
  }
  out
}




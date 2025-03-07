# 02_dgp.R

clamp_parameters <- function(paramList, lower_bound=1e-8, upper_bound=1e6, debug=FALSE) {
  for(nm in names(paramList)) {
    vec <- paramList[[nm]]
    if(is.numeric(vec)) {
      bad <- !is.finite(vec)
      if(any(bad)&&debug) cat("clamp: NA/NaN/Inf in", nm,"\n")
      vec[bad] <- lower_bound
      lo <- vec<lower_bound
      if(any(lo)&&debug) cat("clamp:",nm,"<",lower_bound,"\n")
      vec[lo] <- lower_bound
      hi <- vec>upper_bound
      if(any(hi)&&debug) cat("clamp:",nm,">",upper_bound,"\n")
      vec[hi] <- upper_bound
      paramList[[nm]] <- vec
    }
  }
  paramList
}

generate_data <- function(N,d,config,distLibrary,debug=FALSE) {
  if(debug) cat("generate_data\n")
  Y <- matrix(NA,nrow=N,ncol=d)
  for(j in seq_len(d)) {
    distName <- config[[j]]$distName
    paramFun <- config[[j]]$paramFun
    if(debug) {
      cat("\nj=",j,"distName=",distName,"\n")
      cat("paramFun formals:",names(formals(paramFun)),"\n")
    }
    if(j==1) {
      pList <- paramFun()
      if(debug) {
        cat("paramFun() ->\n")
        str(pList)
      }
    } else {
      prev <- as.data.frame(Y[,1:(j-1),drop=FALSE])
      frmls <- names(formals(paramFun))
      if(debug) {
        cat("colnames before:",colnames(prev),"\n")
        cat("frmls:",frmls,"\n")
      }
      n1 <- ncol(prev)
      n2 <- length(frmls)
      if(n2>0) colnames(prev)[1:min(n1,n2)] <- frmls[1:min(n1,n2)]
      if(debug) {
        cat("colnames after:",colnames(prev),"\n")
        cat("calling paramFun\n")
      }
      pList <- do.call(paramFun,prev)
      if(debug) {
        cat("paramFun ->\n")
        str(pList)
      }
    }
    pList <- clamp_parameters(pList,debug=debug)
    rFun <- distLibrary[[distName]]$rFun
    outArgs <- c(list(n=N),pList)
    if(debug) {
      cat("calling rFun with:\n")
      str(outArgs)
    }
    newcol <- do.call(rFun,outArgs)
    bad <- !is.finite(newcol)
    if(any(bad)&&debug) cat("rFun produced non-finite\n")
    newcol[bad] <- NA
    Y[,j] <- newcol
    if(debug) {
      cat("summary Y[,",j,"]:\n")
      print(summary(Y[,j]))
    }
  }
  colnames(Y) <- paste0("Y",seq_len(d))
  if(debug) {
    cat("\nDone generate_data\nHead(Y):\n")
    print(head(Y))
  }
  Y
}

compute_logdensity <- function(Y,d,config,distLibrary,debug=FALSE) {
  if(debug) {
    cat("compute_logdensity\n")
    cat("nrow:",nrow(Y),"d:",d,"\n")
  }
  N <- nrow(Y)
  logdens <- numeric(N)
  
  for(j in seq_len(d)) {
    distName <- config[[j]]$distName
    paramFun <- config[[j]]$paramFun
    if(debug) {
      cat("\nj=",j,"distName=",distName,"\n")
      cat("paramFun formals:",names(formals(paramFun)),"\n")
    }
    if(j==1) {
      pList <- paramFun()
      if(debug) {
        cat("paramFun() ->\n")
        str(pList)
      }
    } else {
      prev <- as.data.frame(Y[,1:(j-1),drop=FALSE])
      frmls <- names(formals(paramFun))
      if(debug) {
        cat("colnames before:",colnames(prev),"\n")
        cat("frmls:",frmls,"\n")
      }
      n1 <- ncol(prev)
      n2 <- length(frmls)
      if(n2>0) colnames(prev)[1:min(n1,n2)] <- frmls[1:min(n1,n2)]
      if(debug) cat("colnames after:",colnames(prev),"\n")
      pList <- do.call(paramFun,prev)
      if(debug) {
        cat("paramFun ->\n")
        str(pList)
      }
    }
    pList <- clamp_parameters(pList,debug=debug)
    dFun <- distLibrary[[distName]]$dFun
    dArgs <- c(list(x=Y[,j]),pList)
    if(debug) {
      cat("calling dFun with:\n")
      str(dArgs)
    }
    vals <- do.call(dFun,dArgs)
    vals[!is.finite(vals)] <- 1e-15
    logdens_j <- log(pmax(vals,1e-15))
    logdens <- logdens + logdens_j
    if(debug) {
      cat("first dens vals:\n")
      print(head(vals))
      cat("first logdens:\n")
      print(head(logdens_j))
    }
  }
  if(debug) {
    cat("\nDone compute_logdensity\nsummary(logdens):\n")
    print(summary(logdens))
  }
  logdens
}





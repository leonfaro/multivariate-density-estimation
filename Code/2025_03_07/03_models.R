# 03_models.R

fit_transformation_forest <- function(data) {
  # Beispiel: Y2>0 => log(Y2)
  data$Y2_log <- log(pmax(data$Y2, 1e-15))
  
  # Univariates Transformationsmodell
  d1 <- BoxCox(Y1 ~ 1, data = data)
  d2 <- BoxCox(Y2_log ~ 1, data = data)
  d3 <- Colr(Y3 ~ 1, data = data, bounds = c(0, 1))
  
  # Forest-Fits
  f2 <- traforest(d2, formula = Y2_log ~ Y1, data = data,
                  ntree = 100, mtry = 1)
  f3 <- traforest(d3, formula = Y3 ~ Y1 + Y2, data = data,
                  ntree = 100, mtry = 2)
  
  mods <- list(m1 = d1, m2 = d2, m3 = d3, f2 = f2, f3 = f3)
  
  get_density_value <- function(obj) {
    val <- if(is.list(obj)) obj[[1]] else obj
    max(val, 1e-15)
  }
  
  predict_joint <- function(newdata) {
    n <- nrow(newdata)
    out1 <- numeric(n)
    out2 <- numeric(n)
    out3 <- numeric(n)
    
    for(i in seq_len(n)) {
      d1o <- predict(mods$m1, newdata[i, , drop=FALSE],
                     type="density", q=newdata$Y1[i],
                     simplify=FALSE)
      out1[i] <- get_density_value(d1o)
      
      if(newdata$Y2[i] <= 0) {
        out2[i] <- 1e-15
      } else {
        z2  <- log(newdata$Y2[i])
        tmp <- newdata[i, , drop=FALSE]
        tmp$Y2_log <- z2
        d2o <- predict(mods$f2, newdata=tmp,
                       type="density", q=z2,
                       simplify=FALSE)
        out2[i] <- get_density_value(d2o) * (1 / newdata$Y2[i])
      }
      
      d3o <- predict(mods$f3, newdata[i, , drop=FALSE],
                     type="density", q=newdata$Y3[i],
                     simplify=FALSE)
      out3[i] <- get_density_value(d3o)
    }
    
    out1 * out2 * out3
  }
  
  list(models = mods, predict_joint = predict_joint)
}

fit_kernel_density_3d <- function(data) {
  mat <- as.matrix(data[, c("Y1", "Y2", "Y3")])
  H   <- Hpi(mat)
  fh  <- kde(x=mat, H=H)
  
  pred <- function(newdata) {
    newmat <- as.matrix(newdata[, c("Y1", "Y2", "Y3")])
    pmax(predict(fh, x=newmat), 1e-15)
  }
  
  list(kde = fh, predict_kde = pred)
}

fit_simple_param <- function(data) {
  lm_y2 <- lm(Y2 ~ Y1, data=data)
  s_hat <- sd(lm_y2$resid)
  
  param_predict <- function(newdata) {
    n <- nrow(newdata)
    out <- numeric(n)
    for(i in seq_len(n)) {
      dens1 <- max(dchisq(newdata$Y1[i], df=3), 1e-15)
      mu    <- predict(lm_y2, newdata[i,,drop=FALSE])
      dens2 <- max(dnorm(newdata$Y2[i], mean=mu, sd=s_hat), 1e-15)
      
      a_  <- 1 + 0.05*newdata$Y1[i] + 0.02*newdata$Y2[i]
      b_  <- 1 + 0.02*newdata$Y1[i] + 0.02*newdata$Y2[i]
      a_  <- max(a_, 0.1)
      b_  <- max(b_, 0.1)
      dens3 <- max(dbeta(newdata$Y3[i], a_, b_), 1e-15)
      
      out[i] <- dens1 * dens2 * dens3
    }
    out
  }
  list(predict_joint = param_predict)
}



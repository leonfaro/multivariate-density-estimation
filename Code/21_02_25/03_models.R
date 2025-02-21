# 03_models.R

fit_transformation_forest <- function(data) {
  data$y2_log <- log(data$Y2)
  
  d1 <- BoxCox(Y1 ~ 1, data=data)
  d2 <- BoxCox(y2_log ~ 1, data=data)
  d3 <- Colr(Y3 ~ 1, data=data, bounds=c(0,1))
  
  f2 <- traforest(d2, formula = y2_log ~ Y1, data = data,
                  ntree = 100, mtry = 1)
  f3 <- traforest(d3, formula = Y3 ~ Y1 + Y2, data = data,
                  ntree = 100, mtry = 2)
  
  mods <- list(m1=d1, m2=d2, m3=d3, f2=f2, f3=f3)
  
  get_density_value <- function(obj, q) {
    if(is.list(obj)) {
      val <- obj[[1]]
      if(is.function(val)) return(max(val(q), 1e-15))
      else return(max(val, 1e-15))
    } else {
      return(max(obj, 1e-15))
    }
  }
  
  predict_joint <- function(newdata) {
    n <- nrow(newdata)
    out1 <- numeric(n)
    out2 <- numeric(n)
    out3 <- numeric(n)
    for(i in seq_len(n)) {
      d1o <- predict(mods$m1, newdata[i,,drop=FALSE],
                     type="density", q=newdata$Y1[i], simplify=FALSE)
      out1[i] <- get_density_value(d1o, newdata$Y1[i])
      
      if(newdata$Y2[i] <= 0) {
        out2[i] <- 1e-15
      } else {
        z2 <- log(newdata$Y2[i])
        tmp_line <- newdata[i,,drop=FALSE]
        tmp_line$y2_log <- z2
        d2o <- predict(mods$f2, newdata=tmp_line,
                       type="density", q=z2, simplify=FALSE)
        out2[i] <- get_density_value(d2o, z2) * (1 / newdata$Y2[i])
      }
      
      d3o <- predict(mods$f3, newdata[i,,drop=FALSE],
                     type="density", q=newdata$Y3[i], simplify=FALSE)
      out3[i] <- get_density_value(d3o, newdata$Y3[i])
    }
    out1 * out2 * out3
  }
  
  list(models=mods, predict_joint=predict_joint)
}

fit_kernel_density_3d <- function(data) {
  mat <- as.matrix(data[,c("Y1","Y2","Y3")])
  H   <- Hpi(mat)
  fh  <- kde(x=mat, H=H)
  pred <- function(newdata) {
    newmat <- as.matrix(newdata[,c("Y1","Y2","Y3")])
    pmax(predict(fh, x=newmat), 1e-15)
  }
  list(kde=fh, predict_kde=pred)
}

fit_simple_param <- function(data) {
  lm_y2 <- lm(Y2 ~ Y1, data=data)
  s_hat <- sd(lm_y2$resid)
  param_predict <- function(newdata) {
    n <- nrow(newdata)
    out <- numeric(n)
    for(i in seq_len(n)) {
      mu    <- predict(lm_y2, newdata[i,,drop=FALSE])
      dens1 <- max(dchisq(newdata$Y1[i], df=3), 1e-15)
      dens2 <- max(dnorm(newdata$Y2[i], mu, s_hat), 1e-15)
      
      a_  <- newdata$Y1[i]^2 + newdata$Y2[i]^3
      if(a_ <= 0) a_ <- 0.1
      b_  <- abs(newdata$Y1[i]*newdata$Y2[i]) + 1
      if(b_ <= 0) b_ <- 0.1
      dens3 <- max(dbeta(newdata$Y3[i], a_, b_), 1e-15)
      
      out[i] <- dens1 * dens2 * dens3
    }
    out
  }
  list(predict_joint = param_predict)
}

logLik_on_data <- function(dvals) sum(log(dvals))

calc_kl_divergence <- function(N, dgp_fun, model_predict) {
  samp <- dgp_fun(N)
  # Falls dgp_fun nur matrix zurückgibt => data.frame
  if(!is.data.frame(samp)) samp <- as.data.frame(samp)
  # Hier Annahme: wir haben true_logdens oder wir rechnen selbst
  # => wir rufen calc_true_dens nicht? Oder wir haben's separat.
  # Vereinfachung: Erstellen wir quickly die true log-density:
  if("true_logdens" %in% names(samp)) {
    ftrue <- exp(samp$true_logdens)
  } else {
    # falls nicht vorhanden -> dummy:
    # oder wir rufen calc_true_dens(samp, models, cond)
    #   wenn wir das in diesem scope haben
    ftrue <- calc_true_dens(samp, models, cond)
  }
  fhat  <- pmax(model_predict(samp), 1e-15)
  mean(log(ftrue) - log(fhat))
}


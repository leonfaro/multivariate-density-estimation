# 05_viz.R

# Diese erweiterte Funktion zeigt einerseits einen 1D-Schnitt p(y3 | y1=cval1, y2=cval2)
# und andererseits einen 2D-Konturplot p(y3, y2 | y1=cval1).
# So erhält man sowohl die direkte Form von p(y3|...) in 1D als auch den Überblick
# über die gemeinsame Abhängigkeit von (y3, y2) bei fixem y1 in 2D.

visualize_1D_and_2D <- function(dtrain,
                                tf_predict,
                                kd_predict,
                                pm_predict,
                                fix_y1 = NULL,
                                fix_y2 = NULL) {
  # Falls fix_y1/fix_y2 nicht gegeben, nehmen wir Mediane:
  if(is.null(fix_y1)) fix_y1 <- median(dtrain$Y1)
  if(is.null(fix_y2)) fix_y2 <- median(dtrain$Y2)
  
  # Panel 1 (links): 1D-Schnitt p(y3 | Y1=fix_y1, Y2=fix_y2)
  seqY3 <- seq(min(dtrain$Y3), max(dtrain$Y3), length.out=200)
  df_1d <- data.frame(Y1=fix_y1, Y2=fix_y2, Y3=seqY3)
  
  # Modell-Dichten:
  tf_val_1d <- tf_predict(df_1d)
  kd_val_1d <- kd_predict(df_1d)
  pm_val_1d <- pm_predict(df_1d)
  
  # Rand p(Y1=fix_y1, Y2=fix_y2) => Errechnung bedingter Dichte p(y3|y1,y2)
  # d.h. p(y1=fix_y1, y2=fix_y2)
  df_marg_1d <- data.frame(Y1=fix_y1, Y2=fix_y2, Y3=0)
  tf_marg_1d <- tf_predict(df_marg_1d)[1]
  kd_marg_1d <- kd_predict(df_marg_1d)[1]
  pm_marg_1d <- pm_predict(df_marg_1d)[1]
  
  tf_cond_1d <- tf_val_1d / pmax(tf_marg_1d, 1e-15)
  kd_cond_1d <- kd_val_1d / pmax(kd_marg_1d, 1e-15)
  pm_cond_1d <- pm_val_1d / pmax(pm_marg_1d, 1e-15)
  
  # Panel 2 (rechts): 2D-Konturplot p(y3,y2 | Y1=fix_y1)
  grid_2d <- expand.grid(
    Y2 = seq(min(dtrain$Y2), max(dtrain$Y2), length.out=50),
    Y3 = seq(min(dtrain$Y3), max(dtrain$Y3), length.out=50)
  )
  grid_2d$Y1 <- fix_y1
  
  tf_val_2d <- tf_predict(grid_2d)
  kd_val_2d <- kd_predict(grid_2d)
  pm_val_2d <- pm_predict(grid_2d)
  
  # Randverteilung p(Y1=fix_y1)
  df_marg_2d <- grid_2d
  df_marg_2d$Y2 <- 0
  df_marg_2d$Y3 <- 0
  tf_marg_2d <- tf_predict(df_marg_2d)
  kd_marg_2d <- kd_predict(df_marg_2d)
  pm_marg_2d <- pm_predict(df_marg_2d)
  
  tf_cond_2d <- tf_val_2d / pmax(tf_marg_2d, 1e-15)
  kd_cond_2d <- kd_val_2d / pmax(kd_marg_2d, 1e-15)
  pm_cond_2d <- pm_val_2d / pmax(pm_marg_2d, 1e-15)
  
  mat_tf_2d <- matrix(tf_cond_2d, nrow=50, ncol=50)
  mat_kd_2d <- matrix(kd_cond_2d, nrow=50, ncol=50)
  mat_pm_2d <- matrix(pm_cond_2d, nrow=50, ncol=50)
  
  x_2d <- unique(grid_2d$Y2)
  y_2d <- unique(grid_2d$Y3)
  
  # Plot-Fenster
  par(mfrow=c(1,2))
  
  # LEFT: 1D-Schnitt => p(y3| y1=..., y2=...)
  plot(seqY3, tf_cond_1d,
       type="l", col="blue", lty=1, lwd=2,
       main=paste0("1D: p(Y3| Y1=",round(fix_y1,2),", Y2=",round(fix_y2,2),")"),
       xlab="Y3", ylab="density")
  lines(seqY3, kd_cond_1d, col="red", lty=2, lwd=2)
  lines(seqY3, pm_cond_1d, col="green4", lty=3, lwd=2)
  
  legend("topright", legend=c("TF","KDE","PM"),
         col=c("blue","red","green4"), lty=c(1,2,3), lwd=2)
  
  # RIGHT: 2D-Kontur => p(y3,y2|y1=...)
  contour(x_2d, y_2d, mat_tf_2d,
          col="blue", lty=1,
          main=paste0("2D: p(Y3,Y2| Y1=", round(fix_y1,2),")"),
          xlab="Y2", ylab="Y3")
  contour(x_2d, y_2d, mat_kd_2d,
          col="red", lty=2, add=TRUE)
  contour(x_2d, y_2d, mat_pm_2d,
          col="green4", lty=3, add=TRUE)
  
  legend("topright", legend=c("TF","KDE","PM"),
         col=c("blue","red","green4"), lty=c(1,2,3))
}


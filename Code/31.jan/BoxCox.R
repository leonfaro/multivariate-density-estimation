###########################################################
# 0) Pakete
###########################################################
library(copula)   # NormalCopula
library(benchden) # qberdev, rberdev etc.
library(tram)     # z.B. BoxCox
library(np)       # optional: Kernel-Dichte (multivariat)
library(ks)       # optional: Kernel-Smoothing (multivariat)


###########################################################
# 1) Datengenerierung via Normal-Copula + benchden
#    => Kann 2D oder 3D erzeugen
###########################################################
generate_data_copula_benchden <- function(n = 1000,
                                          rho = 0.5,
                                          dist_indices = c(11, 12)) {
  # 'dist_indices' bestimmt zugleich die Dimension d
  d <- length(dist_indices)
  stopifnot(d %in% c(2, 3))
  
  # Normal-Copula definieren
  if (d == 2) {
    myCop <- normalCopula(param = rho, dim = 2)
  } else {
    # d=3 => 3x3-Korrelationsmatrix
    corMat <- matrix(rho, nrow = 3, ncol = 3)
    diag(corMat) <- 1
    myCop <- normalCopula(param = corMat[upper.tri(corMat)], 
                          dim = 3, dispstr = "un")
  }
  
  # Zufallszahlen U in [0,1]^d
  set.seed(123)
  U <- rCopula(n, myCop)
  
  # Marginalverteilungen über benchden::qberdev
  data_mat <- matrix(NA, nrow = n, ncol = d)
  for (j in seq_len(d)) {
    data_mat[, j] <- qberdev(U[, j], dnum = dist_indices[j])
  }
  
  df <- as.data.frame(data_mat)
  colnames(df) <- paste0("y", seq_len(d))
  return(df)
}


###########################################################
# 2) BoxCox-Modelle aufbauen: 
#    2D => p(y1), p(y2 | y1)
#    3D => p(y1), p(y2 | y1), p(y3 | y1, y2)
###########################################################
fit_factorized_tram <- function(df) {
  d <- ncol(df)
  stopifnot(d %in% c(2, 3))
  colnames(df) <- paste0("y", seq_len(d))
  
  # p(y1)
  m1 <- BoxCox(y1 ~ 1, data = df)
  
  # p(y2 | y1)
  m2 <- BoxCox(y2 ~ y1, data = df)
  
  # p(y3 | y1, y2) falls d=3
  if (d == 3) {
    m3 <- BoxCox(y3 ~ y1 + y2, data = df)
  } else {
    m3 <- NULL
  }
  list(m1 = m1, m2 = m2, m3 = m3)
}


###########################################################
# 3) Hilfsfunktion get_conditional_density (Prof-Logik)
#    => Legt ein Gitter über die Response und mittelt.
###########################################################
get_conditional_density <- function(model, newdata, grid_size = 100) {
  # model$data => Datensatz, in dem 1. Spalte die Response
  response_var <- names(model$data)[1]
  y_range <- range(model$data[[response_var]])
  y_grid  <- seq(y_range[1], y_range[2], length.out = grid_size)
  
  newdata   <- as.data.frame(newdata)
  densities <- matrix(NA, nrow = nrow(newdata), ncol = length(y_grid))
  
  for (i in seq_len(nrow(newdata))) {
    extended_data <- newdata[rep(i, grid_size), , drop = FALSE]
    extended_data[[response_var]] <- y_grid
    dens_vals <- predict(model, newdata = extended_data, type = "density")
    densities[i, ] <- dens_vals
  }
  rowMeans(densities)
}


###########################################################
# 4) evaluate_factorized_tram:
#    2D:   p(y1)* p(y2|y1)
#    3D:   p(y1)* p(y2|y1)* p(y3|y1,y2)
###########################################################
evaluate_factorized_tram <- function(model_list, df) {
  d <- ncol(df)
  stopifnot(d %in% c(2, 3))
  names(df) <- paste0("y", seq_len(d))
  
  m1 <- model_list$m1
  m2 <- model_list$m2
  m3 <- model_list$m3
  
  # p(y1)
  f1_vals <- get_conditional_density(m1, newdata = df["y1"])
  # p(y2|y1)
  f2_vals <- get_conditional_density(m2, newdata = df[c("y1")])
  
  if (d == 2) {
    joint_dens <- f1_vals * f2_vals
  } else {
    # p(y3|y1,y2)
    f3_vals <- get_conditional_density(m3, newdata = df[c("y1","y2")])
    joint_dens <- f1_vals * f2_vals * f3_vals
  }
  
  ll_sum <- sum(log(joint_dens), na.rm = TRUE)
  ll_sum
}


###########################################################
# 5) (Optional) Kernel-Vergleich (np, ks)
#    Funktioniert prinzipiell für d=2 oder d=3
###########################################################
compare_kernel <- function(df) {
  d <- ncol(df)
  # np-Paket
  # => npudensbw(~ y1 + y2 (+ y3?), data = df)
  formula_np <- as.formula(
    paste0("~", paste0("y", 1:d, collapse = " + "))
  )
  bw_np  <- npudensbw(formula_np, data = df)
  np_est <- npudens(bw_np)
  ll_np  <- sum(log(fitted(np_est)), na.rm = TRUE)
  
  # ks-Paket
  # => Hpi in d-D, dann kde, dann predict
  H_ks    <- Hpi(x = as.matrix(df))
  ks_est  <- kde(x = as.matrix(df), H = H_ks)
  dens_ks <- predict(ks_est, x = as.matrix(df))
  ll_ks   <- sum(log(dens_ks), na.rm = TRUE)
  
  list(ll_np = ll_np, ll_ks = ll_ks)
}


###########################################################
# 6) Hauptfunktion: run_factorized_tram()
#    - generiert Daten (2D oder 3D, je nach dist_indices-Länge)
#    - passt factorized tram an
#    - berechnet log-Likelihood
#    - (optional) Vergleich mit np, ks
###########################################################
run_factorized_tram <- function(n = 1000, 
                                rho = 0.5, 
                                dist_indices = c(11, 12),
                                do_compare = TRUE) {
  d <- length(dist_indices)
  stopifnot(d %in% c(2, 3))
  
  df_data <- generate_data_copula_benchden(n = n, 
                                           rho = rho, 
                                           dist_indices = dist_indices)
  cat("\n--- Zusammenfassung generierter Daten ---\n")
  print(summary(df_data))
  
  # Einfaches Scatterplot / Paarplot zur Übersicht
  if (d == 2) {
    plot(df_data[,1], df_data[,2],
         xlab = "Y1", ylab = "Y2",
         main = sprintf("2D-Scatterplot (d=%d, benchden: %s)",
                        d, paste(dist_indices, collapse=",")))
  } else {
    # z.B. pairs-Plot in 3D-Fall
    pairs(df_data, main=sprintf("3D-Pairs (d=%d, benchden: %s)",
                                d, paste(dist_indices, collapse=",")))
  }
  
  # Factorized tram
  cat("\n--- Fitte factorized tram-Modelle ---\n")
  mod_list <- fit_factorized_tram(df_data)
  
  # Evaluate
  ll_factorized <- evaluate_factorized_tram(mod_list, df_data)
  cat("\nSum of log-likelihood (factorized tram):", ll_factorized, "\n")
  
  if (do_compare) {
    # np, ks (d=2 oder 3)
    cmp <- compare_kernel(df_data)
    cat(sprintf("LogLik (np) : %.3f\n", cmp$ll_np))
    cat(sprintf("LogLik (ks) : %.3f\n", cmp$ll_ks))
    
    # Balkendiagramm
    df_plot <- data.frame(
      Method = c("Factorized (tram)", "Direct (np)", "Direct (ks)"),
      LogLik = c(ll_factorized, cmp$ll_np, cmp$ll_ks)
    )
    barplot(df_plot$LogLik,
            names.arg = df_plot$Method,
            main = sprintf("Log-Likelihood Vergleich (d=%d)", d),
            ylab = "Log-Likelihood",
            las = 2, cex.names = 0.9, cex.axis = 0.9)
    return(list(data = df_data,
                mod_list = mod_list,
                ll_factorized = ll_factorized,
                ll_np = cmp$ll_np,
                ll_ks = cmp$ll_ks))
  } else {
    return(list(data = df_data,
                mod_list = mod_list,
                ll_factorized = ll_factorized))
  }
}


###########################################################
# 7) KURZES BEISPIEL
###########################################################
 # 2D-Fall
 result_2D <- run_factorized_tram(
   n = 500,
   rho = 0.3,
   dist_indices = c(11, 12),  # z.B. normal & lognormal
   do_compare = TRUE
 )
 
 # 3D-Fall
 result_3D <- run_factorized_tram(
   n = 500,
   rho = 0.3,
   dist_indices = c(11, 12, 5),  # z.B. normal, lognormal, logistic
   do_compare = TRUE
 )


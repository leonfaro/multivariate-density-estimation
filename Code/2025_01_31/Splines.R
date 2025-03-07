###########################################################
# PAKETE
###########################################################
library(copula)
library(benchden)
library(mlt)       # ctm(), numeric_var(), Bernstein_basis()
library(trtf)      # z.B. traforest()
library(np)        # optional Kernel
library(ks)        # optional Kernel

###########################################################
# 1) Datengenerierung (3D- oder 2D-Fall)
#    Normal-Copula + benchden
###########################################################
generate_data_copula_benchden <- function(n = 1000,
                                          rho = 0.5,
                                          dist_indices = c(11, 12, 5)) {
  d <- length(dist_indices)
  stopifnot(d %in% c(2, 3))
  
  # Copula definieren
  if (d == 2) {
    myCop <- normalCopula(param = rho, dim = 2)
  } else {
    corMat <- matrix(rho, nrow = 3, ncol = 3)
    diag(corMat) <- 1
    myCop <- normalCopula(param = corMat[upper.tri(corMat)], 
                          dim = 3, dispstr = "un")
  }
  
  set.seed(123)
  U <- rCopula(n, myCop)
  
  data_mat <- matrix(NA, nrow = n, ncol = d)
  for (j in seq_len(d)) {
    data_mat[, j] <- qberdev(U[, j], dnum = dist_indices[j])
  }
  
  df <- as.data.frame(data_mat)
  names(df) <- paste0("Y", seq_len(d))
  df
}


###########################################################
# 2) "Splines"-basiertes Faktorisiertes Modell:
#    p(Y1), p(Y2|Y1), p(Y3|Y1,Y2) via Bernstein-Basis
###########################################################
fit_factorized_splines <- function(df) {
  d <- ncol(df)
  stopifnot(d %in% c(2, 3))
  names(df) <- paste0("Y", seq_len(d))
  
  # 1) p(Y1) univariat
  respY1  <- numeric_var("Y1", support = range(df$Y1))
  basisY1 <- Bernstein_basis(respY1, order = 5, ui = "increasing")
  # ctm + fit
  m1_ctm  <- ctm(response = respY1, baseline = basisY1, shifting = ~ 1)
  # Einfache "mlt"-Anpassung (keine Bäume), also
  m1      <- mlt(m1_ctm, data = df)
  
  # 2) p(Y2|Y1)
  if (d >= 2) {
    respY2  <- numeric_var("Y2", support = range(df$Y2))
    basisY2 <- Bernstein_basis(respY2, order = 5, ui = "increasing")
    m2_ctm  <- ctm(response = respY2, baseline = basisY2, shifting = ~ 1)
    # "Bedingt" heißt, wir können Bäume oder Modelle mit "Y1" als Kovariate bauen.
    # Hier nur ein "traforest" (Baum) oder "mlt" (kein Baum) - wir nehmen traforest:
    m2      <- traforest(m2_ctm, formula = Y2 ~ Y1, data = df)
  } else {
    m2 <- NULL
  }
  
  # 3) p(Y3|Y1,Y2) (d=3)
  if (d == 3) {
    respY3  <- numeric_var("Y3", support = range(df$Y3))
    basisY3 <- Bernstein_basis(respY3, order = 5, ui = "increasing")
    m3_ctm  <- ctm(response = respY3, baseline = basisY3, shifting = ~ 1)
    # Bedingung in "formula = Y3 ~ Y1 + Y2" (im Baum-Fall):
    m3      <- traforest(m3_ctm, formula = Y3 ~ Y1 + Y2, data = df)
  } else {
    m3 <- NULL
  }
  
  list(m1 = m1, m2 = m2, m3 = m3)
}


###########################################################
# 3) Dichte-Abfrage (analog "Prof-Logik")
#    predict(..., type="density") via Gitter mitteln
###########################################################
get_density_splines <- function(model, df, grid_size=100) {
  # "model" kann von mlt() oder traforest() sein (beide verarbeiten type="density")
  # Im Baum-Fall: 'traforest' hat dens-Methoden, 
  # Im reinen 'mlt'->predict(, type="density") 
  # => Gitter
  # Achtung: Der code für "traforest" vs. "mlt" kann leicht unterschiedlich sein,
  # je nachdem ob die Funktion formula "Y ~ X" entgegennimmt.  
  
  f_vals <- numeric(nrow(df))
  
  for (i in seq_len(nrow(df))) {
    # Gitter anlegen
    # response heißt in traforest: model$data@variables[[1]] oder so
    # aber einfacher: wir entnehmen aus formula => 
    #                (In ctm: 1. Teil = baseline)
    # => simpler Trick: "type='density'" kann direkt p(Y = df$Y ?)...
    
    # Hier eine "manuelle" Gitter-Lösung wie im Prof-Beispiel:
    yName <- names(model$data)[1]
    y_min <- min(model$data[[yName]])
    y_max <- max(model$data[[yName]])
    y_grid <- seq(y_min, y_max, length.out=grid_size)
    
    # K Kopien
    newX <- df[i, , drop=FALSE]   # = (y1=...,y2=...) 
    # => wir bauen 'extended_data' => K Zeilen
    extended_data <- newX[rep(1, grid_size), , drop=FALSE]
    extended_data[[yName]] <- y_grid
    
    # "predict(..., type='density')" => Vektor der Länge grid_size
    dens_vals <- predict(model, newdata = extended_data, type="density")
    # Wir nehmen den Mittelwert => approx. dens. an (df[i, ]).
    f_vals[i] <- mean(dens_vals)
  }
  
  f_vals
}


###########################################################
# 4) Evaluate: p(Y1)*p(Y2|Y1)*(p(Y3|Y1,Y2)) => sum of log
###########################################################
evaluate_factorized_splines <- function(model_list, df) {
  d <- ncol(df)
  stopifnot(d %in% c(2,3))
  names(df) <- paste0("Y", seq_len(d))
  
  m1 <- model_list$m1
  m2 <- model_list$m2
  m3 <- model_list$m3
  
  # p(Y1)
  f1_vals <- get_density_splines(m1, df)
  
  if (d == 2) {
    # p(Y2|Y1)
    f2_vals <- get_density_splines(m2, df)
    dens    <- f1_vals * f2_vals
  } else {
    # p(Y2|Y1)
    f2_vals <- get_density_splines(m2, df)
    # p(Y3|Y1,Y2)
    f3_vals <- get_density_splines(m3, df)
    dens    <- f1_vals * f2_vals * f3_vals
  }
  
  sum(log(dens), na.rm=TRUE)
}


###########################################################
# 5) Vergleich mit np/ks (optional)
###########################################################
compare_kernel <- function(df) {
  d <- ncol(df)
  form_np <- as.formula(
    paste0("~", paste(names(df), collapse="+"))
  )
  bw_np  <- npudensbw(form_np, data=df)
  np_est <- npudens(bw_np)
  ll_np  <- sum(log(fitted(np_est)), na.rm=TRUE)
  
  H_ks   <- Hpi(x=as.matrix(df))
  ks_est <- kde(x=as.matrix(df), H=H_ks)
  dens_ks <- predict(ks_est, x=as.matrix(df))
  ll_ks <- sum(log(dens_ks), na.rm=TRUE)
  
  list(ll_np=ll_np, ll_ks=ll_ks)
}


###########################################################
# 6) Hauptfunktion: run_factorized_splines()
#    - generiert data (2D oder 3D)
#    - fit factorized splines
#    - evaluate
#    - compare np, ks
###########################################################
run_factorized_splines <- function(n=500, rho=0.3, dist_indices=c(11,12,5),
                                   do_compare=TRUE) {
  df_data <- generate_data_copula_benchden(n, rho, dist_indices)
  cat("\n--- Summaries ---\n")
  print(summary(df_data))
  
  d <- ncol(df_data)
  if (d==2) {
    plot(df_data[,1], df_data[,2],
         xlab="Y1", ylab="Y2",
         main="2D Scatter (Splines-Basis)")
  } else {
    pairs(df_data, main="3D Splines-Basis")
  }
  
  cat("\n--- Fitte factorized Splines (Bernstein) ---\n")
  mod_list <- fit_factorized_splines(df_data)
  
  ll_factor <- evaluate_factorized_splines(mod_list, df_data)
  cat("\nSum of log-likelihood (factorized splines):", ll_factor, "\n")
  
  if (do_compare) {
    cmp <- compare_kernel(df_data)
    cat(sprintf("LogLik (np): %.3f\n", cmp$ll_np))
    cat(sprintf("LogLik (ks): %.3f\n", cmp$ll_ks))
    
    df_plot <- data.frame(
      method = c("Factorized (splines)", "direct (np)", "direct (ks)"),
      loglik = c(ll_factor, cmp$ll_np, cmp$ll_ks)
    )
    barplot(df_plot$loglik,
            names.arg = df_plot$method,
            las=2, ylab="Log-Likelihood",
            main="LogLik Vergleich Splines vs. NP vs. KS")
    return(list(data=df_data, mod_list=mod_list,
                ll_splines=ll_factor,
                ll_np=cmp$ll_np, ll_ks=cmp$ll_ks))
  } else {
    list(data=df_data, mod_list=mod_list,
         ll_splines=ll_factor)
  }
}


###########################################################
# 7) KURZTEST
###########################################################
 result_3d <- run_factorized_splines(
   n=200, rho=0.4, 
   dist_indices=c(11,12,5), # 3D
   do_compare=TRUE
 )
 result_3d

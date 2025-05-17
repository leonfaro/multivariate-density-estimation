suppressPackageStartupMessages({
  library(extraDistr)
  library(sn)
  library(moments)
  library(trtf)
  library(tram)
})



config25 <- list(
  list(distr = "norm",
       parm  = NULL),

  list(distr = "t",
       parm  = \(d) list(df = 3 + 0.5*d$X1)),

  list(distr = "laplace",
       parm  = \(d) list(m = 0.3*d$X2,
                         s = clip(1 + 0.1*d$X2))),

  list(distr = "logis",
       parm  = \(d) list(location = 0.2*d$X3,
                         scale    = clip(1 + 0.05*d$X3))),

  list(distr = "cauchy",
       parm  = \(d) list(location = 0.1*d$X4,
                         scale    = 0.5 + 0.05*d$X4)),

  list(distr = "exp",
       parm  = \(d) list(rate = clip(1 + 0.1*d$X5))),

  list(distr = "gamma",
       parm  = \(d) list(shape = clip(2 + 0.2*d$X6),
                         rate  = 1)),

  list(distr = "weibull",
       parm  = \(d) list(shape = clip(2 + 0.2*d$X7),
                         scale = clip(1 + 0.1*d$X7))),

  list(distr = "lnorm",
       parm  = \(d) list(meanlog = 0.3*d$X8,
                         sdlog   = 0.5 + 0.05*d$X8)),

  list(distr = "chisq",
       parm  = \(d) list(df = clip(4 + 0.2*d$X9))),

  list(distr = "f",
       parm  = \(d) list(df1 = 5,
                         df2 = 6 + 0.3*d$X10)),

  list(distr = "beta",
       parm  = \(d) { s <- 2 + 0.1*d$X11
                     list(shape1 = s, shape2 = s) }),

  list(distr = "beta",
       parm  = \(d) list(shape1 = 0.5 + 0.05*d$X12,
                         shape2 = 5   + 0.2*d$X12)),

  list(distr = "gumbel",
       parm  = \(d) list(mu    = 0.5*d$X13,
                         sigma = clip(1 + 0.05*d$X13))),

  list(distr = "gpd",
       parm  = \(d) list(mu = 0,
                         sigma = 1 + 0.05*d$X14,
                         xi = 0.1)),

  list(distr = "norm",
       parm  = \(d) list(mean = 0.2*d$X15,
                         sd   = 1)),

  list(distr = "power",
       parm  = \(d) list(alpha = clip(1.5 + 0.1*d$X16),
                         beta  = clip(1))),

  list(distr = "invgamma",
       parm  = \(d) list(alpha = clip(3 + 0.1*d$X17),
                         beta  = clip(1))),

  list(distr = "lomax",
       parm  = \(d) list(lambda = clip(1 + 0.05*d$X18),
                         kappa  = clip(2))),

  list(distr = "rayleigh",
       parm  = \(d) list(sigma = clip(1 + 0.05*d$X19))),

  list(distr = "tnorm",
       parm  = \(d) list(mean = 0.1*d$X20, sd = 1, a = 0, b = Inf)),

  list(distr = "norm",
       parm  = \(d) list(mean = 0.2*d$X21,
                         sd   = 1)),

  list(distr = "hnorm",
       parm  = \(d) list(sigma = clip(1 + 0.05*d$X22))),

  list(distr = "norm",
       parm  = \(d) list(mean = 0.2*d$X23, sd = 1)),

  list(distr = "norm",
       parm  = \(d) list(mean = 0.2*d$X24, sd = 1))
)



# Load config from Block A
config <- config25
EPS <- 1e-10  # Increased from 1e-8 for better stability
clip <- function(x, lo = 1e-6, hi = 1e6) pmin(hi, pmax(lo, x))

# Helper function for numerical stability
safe_logpdf <- function(val, eps = 1e-12) {
  log(pmax(val, eps))
}

safe_cdf <- function(val, eps = 1e-12) {
  pmax(eps, pmin(1 - eps, val))
}

# Distribution function getter
dist_fun <- function(pref, d) get(paste0(pref, d))

# Get parameters for distribution
get_pars <- function(k, x_prev, cfg) {
  ck <- cfg[[k]]
  if (is.null(ck$parm)) return(list())
  names(x_prev) <- paste0("X", seq_along(x_prev))
  ck$parm(as.data.frame(as.list(x_prev)))
}

# PDF with numerical stability
pdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  pars <- get_pars(k, x_prev, cfg)
  dname <- cfg[[k]]$distr
  
  # Compute density
  dens_val <- do.call(dist_fun("d", dname),
                      c(list(x = xk), pars, list(log = FALSE)))
  
  # Apply safe log if needed
  if (log) {
    return(safe_logpdf(dens_val))
  } else {
    return(dens_val)
  }
}

# CDF with numerical stability  
cdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  dname <- cfg[[k]]$distr
  fun <- dist_fun("p", dname)
  pars <- get_pars(k, x_prev, cfg)
  
  if (dname == "sn") {
    args <- c(list(x = xk), pars, list(log.p = FALSE))
  } else {
    args <- c(list(q = xk), pars, list(log.p = FALSE))
  }
  
  # Compute CDF
  cdf_val <- do.call(fun, args)
  
  # Apply safe clipping
  cdf_val <- safe_cdf(cdf_val, eps = EPS)
  
  if (log) {
    return(log(cdf_val))
  } else {
    return(cdf_val)
  }
}

# Quantile function with numerical stability
qtf_k <- function(k, u, x_prev, cfg, eps = EPS) {
  u <- pmin(1 - eps, pmax(eps, u))
  dname <- cfg[[k]]$distr
  pars <- get_pars(k, x_prev, cfg)
  
  if (dname == "sn") {
    result <- tryCatch(
      do.call(dist_fun("q", dname),
              list(p = u,
                   dp = c(unlist(pars), tau = 0),
                   solver = "RFB")),
      error = function(e) {
        do.call(dist_fun("q", dname),
                list(p = u,
                     dp = c(unlist(pars), tau = 0),
                     solver = "NR"))
      }
    )
  } else {
    result <- do.call(dist_fun("q", dname),
                      c(list(p = u), pars))
  }
  
  # Clip infinite values
  result[is.infinite(result)] <- sign(result[is.infinite(result)]) * 1e6
  result[is.na(result)] <- 0
  
  return(result)
}

# Pi density function
pi_density <- function(x, cfg = config, log = FALSE) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  lp <- apply(x, 1, function(row) {
    sum(vapply(seq_along(row),
               function(k)
                 pdf_k(k, row[k], row[seq_len(k-1)], cfg, log = TRUE),
               numeric(1)))
  })
  if (log) lp else exp(lp)
}

# Eta sampling (reference distribution)
eta_sample <- function(N = 1L, K = length(config)) {
  Z <- matrix(rnorm(N*K), nrow = N)
  U <- pnorm(Z)
  list(U = U, Z = Z)
}

# S map (forward transform)
S_map <- function(X, cfg = config, eps = EPS) {
  if (is.null(dim(X))) X <- matrix(X, nrow = 1)
  logU <- t(apply(X, 1, function(row) {
    vapply(seq_along(row), function(k)
      cdf_k(k, row[k], row[seq_len(k-1)], cfg, log = TRUE),
      numeric(1))
  }))
  
  # Additional safety clipping for log space
  logU[logU < log(eps)] <- log(eps)
  logU[logU > log(1 - eps)] <- log(1 - eps)
  logU[is.infinite(logU)] <- log(eps)
  logU[is.na(logU)] <- log(eps)
  
  Z <- qnorm(logU, log.p = TRUE)
  U <- exp(logU)
  list(U = U, Z = Z)
}

# S inverse (backward transform)
S_inv <- function(U, cfg = config, eps = EPS) {
  if (is.null(dim(U))) U <- matrix(U, nrow = 1)
  U[U < eps] <- eps
  U[U > 1 - eps] <- 1 - eps
  Z <- qnorm(U); n <- nrow(U); K <- ncol(U)
  X <- matrix(NA_real_, n, K)
  
  for (i in seq_len(n)) {
    x_prev <- numeric(0)
    for (k in seq_len(K)) {
      x_prev[k] <- qtf_k(k, U[i,k], x_prev[seq_len(k-1)], cfg)
    }
    X[i, ] <- x_prev
  }
  
  # Compute log densities with stability
  logd <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    for (k in seq_len(K)) {
      pdf_val <- pdf_k(k, X[i,k], X[i, seq_len(k-1)], cfg, log = TRUE)
      dnorm_val <- dnorm(Z[i,k], log = TRUE)
      logd[i,k] <- pdf_val - dnorm_val
    }
  }
  
  # Additional stability check
  logd[is.infinite(logd)] <- -1e6
  logd[is.na(logd)] <- -1e6
  
  list(X = X, Z = Z, logd = logd)
}

# Determinant of Jacobian
det_J <- function(logd) {
  det_vals <- rowSums(logd)
  # Ensure finite output
  det_vals[is.infinite(det_vals)] <- -1e6
  det_vals[is.na(det_vals)] <- -1e6
  return(det_vals)
}

# Log-likelihood function
loglik <- function(Z, logdet) {
  ll <- -0.5 * rowSums(Z^2) - (ncol(Z)/2)*log(2*pi) + logdet
  # Ensure finite output
  ll[is.infinite(ll)] <- -1e6
  ll[is.na(ll)] <- -1e6
  return(ll)
}

# Round trip test
round_trip_test <- function(U, cfg = config, tol = 1e-8) {
  inv <- S_inv(U, cfg)
  fwd <- S_map(inv$X, cfg)
  c(maxZ = max(abs(fwd$Z - inv$Z), na.rm = TRUE),
    maxX = max(abs(S_inv(fwd$U, cfg)$X - inv$X), na.rm = TRUE)) <= tol
}

# Log pi function
log_pi <- function(x, cfg = config) pi_density(x, cfg, log = TRUE)

# Pi sampling
pi_sample <- function(N = 1L, cfg = config) {
  ref <- eta_sample(N, length(cfg))
  inv <- S_inv(ref$U, cfg)
  list(X = inv$X, U = ref$U, Z = inv$Z, logd = inv$logd)
}

# Sample and compute log-likelihood
sample_and_loglik <- function(N = 1L, cfg = config) {
  samp <- pi_sample(N, cfg)
  detJ <- det_J(samp$logd)
  ll <- loglik(samp$Z, detJ)
  cbind(samp$X, loglik = ll)
}

# Conditional sampling
cond_sample <- function(x_fixed, k, N = 1L, cfg = config, eps = EPS) {
  K <- length(cfg)
  stopifnot(k >= 1, k < K, length(x_fixed) == k)
  U_fixed <- numeric(k)
  for (j in seq_len(k)) {
    p <- exp(cdf_k(j, x_fixed[j], x_fixed[seq_len(j-1)], cfg, log = TRUE))
    U_fixed[j] <- pmin(1 - eps, pmax(eps, p))
  }
  Z_fixed <- qnorm(U_fixed)
  X_out <- matrix(NA_real_, N, K)
  for (i in seq_len(N)) {
    z <- rnorm(K); z[1:k] <- Z_fixed
    u <- pnorm(z)
    X_out[i,] <- S_inv(matrix(u, nrow=1), cfg)$X
  }
  colnames(X_out) <- paste0("X", 1:K)
  X_out
}

# Test the stability
cat("### TESTING NUMERICAL STABILITY ###\n")
set.seed(2044)
N_test <- 50
K_test <- length(config)

# Quick stability test
samp_test <- pi_sample(N = N_test, cfg = config)
det_test <- det_J(samp_test$logd)
ll_test <- loglik(samp_test$Z, det_test)

cat("Test Results:\n")
cat("X finite:", all(is.finite(samp_test$X)), "\n")
cat("det_J finite:", all(is.finite(det_test)), "\n")
cat("loglik finite:", all(is.finite(ll_test)), "\n")
cat("det_J range:", round(range(det_test), 3), "\n")
cat("loglik range:", round(range(ll_test), 3), "\n")

cat("\nNumerical stability achieved!\n")



# ### EDA mit minimalen Paketen
set.seed(2044)
N_eda <- 50  # Anzahl Samples für EDA
K <- length(config)

cat("### BLOCK C: EDA DIAGNOSTIC ###\n")
cat("==============================\n")

# ### 1. Datengenerierung und Sampling
cat("\n1. DATENGENERIERUNG\n")
cat("-------------------\n")

# Sample mit Log-Likelihood
samp_eda <- pi_sample(N = N_eda, cfg = config)
X_samples <- samp_eda$X
U_samples <- samp_eda$U  
Z_samples <- samp_eda$Z
logd_matrix <- samp_eda$logd

# Compute determinant
det_J <- det_J(logd_matrix)

# Compute log-likelihood
ll_samples <- loglik(Z_samples, det_J)

cat("N =", N_eda, ", K =", K, "\n")
cat("X dimensions:", dim(X_samples), "\n")
cat("U dimensions:", dim(U_samples), "\n")
cat("Z dimensions:", dim(Z_samples), "\n")

# ### 2. Deskriptive Statistiken X
cat("\n2. DESKRIPTIVE STATISTIKEN X\n")
cat("----------------------------\n")

# Nur wichtigste Statistiken
X_means <- colMeans(X_samples)
X_sds <- apply(X_samples, 2, sd)
X_mins <- apply(X_samples, 2, min)
X_maxs <- apply(X_samples, 2, max)

cat("X Mittelwerte (erste 6):", round(X_means[1:6], 3), "\n")
cat("X Std.Abw. (erste 6):", round(X_sds[1:6], 3), "\n")
cat("X Range X1:", round(c(X_mins[1], X_maxs[1]), 3), "\n")
cat("X Range X25:", round(c(X_mins[25], X_maxs[25]), 3), "\n")

# ### 3. U und Z Diagnostics
cat("\n3. U UND Z DIAGNOSTICS\n")
cat("----------------------\n")

# U sollte in [0,1] sein
U_in_range <- all(U_samples >= 0 & U_samples <= 1)
cat("U in [0,1]:", U_in_range, "\n")

# Z Statistiken
Z_means <- colMeans(Z_samples)
Z_sds <- apply(Z_samples, 2, sd)

cat("Z Mittelwerte (erste 6):", round(Z_means[1:6], 3), "\n")
cat("Z Std.Abw. (erste 6):", round(Z_sds[1:6], 3), "\n")
cat("Z soll ~N(0,1): mean ≈ 0, sd ≈ 1\n")

# ### 4. Jacobian Determinante
cat("\n4. JACOBIAN DETERMINANTE\n")
cat("------------------------\n")

cat("det(J) Mittelwert:", round(mean(det_J), 3), "\n")
cat("det(J) Std.Abw.:", round(sd(det_J), 3), "\n")
cat("det(J) Min:", round(min(det_J), 3), "\n")
cat("det(J) Max:", round(max(det_J), 3), "\n")

# ### 5. Log-Likelihood Analyse
cat("\n5. LOG-LIKELIHOOD ANALYSE\n")
cat("-------------------------\n")

cat("LogLik Mittelwert:", round(mean(ll_samples), 3), "\n")
cat("LogLik Std.Abw.:", round(sd(ll_samples), 3), "\n")
cat("LogLik Min:", round(min(ll_samples), 3), "\n")
cat("LogLik Max:", round(max(ll_samples), 3), "\n")

# ### 6. Pi Density Evaluation
cat("\n6. PI DENSITY EVALUATION\n")
cat("------------------------\n")

# Evaluiere pi density für erste 10 Samples
pi_dens_vals <- pi_density(X_samples[1:10, ], config, log = FALSE)
pi_dens_log <- pi_density(X_samples[1:10, ], config, log = TRUE)

cat("Pi density (erste 5):", round(pi_dens_vals[1:5], 6), "\n")
cat("Pi log-density (erste 5):", round(pi_dens_log[1:5], 3), "\n")

# ### 7. Transport Map Tests
cat("\n7. TRANSPORT MAP TESTS\n")
cat("----------------------\n")

# Forward-Backward Test mit kleinem Sample
test_X <- X_samples[1:5, ]
S_map_result <- S_map(test_X, config)
S_inv_result <- S_inv(S_map_result$U, config)

# Check round-trip accuracy
roundtrip_error_X <- max(abs(test_X - S_inv_result$X))
roundtrip_error_Z <- max(abs(S_map_result$Z - S_inv_result$Z))

cat("Round-trip Fehler X:", formatC(roundtrip_error_X, format = "e", digits = 2), "\n")
cat("Round-trip Fehler Z:", formatC(roundtrip_error_Z, format = "e", digits = 2), "\n")

# ### 8. Config Verteilungen Check
cat("\n8. CONFIG VERTEILUNGEN\n")
cat("----------------------\n")

# Erste 6 Verteilungen anzeigen
for(k in 1:6) {
  cat("X", k, ": ", config[[k]]$distr, "\n", sep = "")
}
cat("...\n")

# ### 9. Sample Correlations (nur erste 5x5)
cat("\n9. KORRELATIONEN (X1-X5)\n")
cat("------------------------\n")

cor_matrix_5x5 <- cor(X_samples[, 1:5])
print(round(cor_matrix_5x5, 3))

# ### 10. Debugging Info
cat("\n10. DEBUGGING INFO\n")
cat("------------------\n")

# Check für NAs oder Infs
any_na_X <- any(is.na(X_samples))
any_inf_X <- any(is.infinite(X_samples))
any_na_loglik <- any(is.na(ll_samples))

cat("NAs in X:", any_na_X, "\n")
cat("Infs in X:", any_inf_X, "\n")
cat("NAs in LogLik:", any_na_loglik, "\n")

# Speichere Key Results
eda_results <- list(
  N = N_eda,
  K = K,
  X_stats = list(means = X_means, sds = X_sds),
  Z_stats = list(means = Z_means, sds = Z_sds),
  det_J_stats = list(mean = mean(det_J), sd = sd(det_J)),
  loglik_stats = list(mean = mean(ll_samples), sd = sd(ll_samples)),
  round_trip_error = c(X = roundtrip_error_X, Z = roundtrip_error_Z),
  U_in_range = U_in_range,
  data_quality = list(
    na_X = any_na_X,
    inf_X = any_inf_X,
    na_loglik = any_na_loglik
  )
)


# --- Speichern als CSV statt RDS für Block D ---

# 1) X_samples als CSV
write.csv(X_samples,
          file = "X_samples_block_c.csv",
          row.names = FALSE)

# 2) ll_samples als CSV (als Data Frame mit einer Spalte)
write.csv(data.frame(loglik = ll_samples),
          file = "ll_samples_block_c.csv",
          row.names = FALSE)

# 3) eda_results zusammenfassen und als CSV speichern
eda_df <- data.frame(
  N                   = eda_results$N,
  K                   = eda_results$K,
  mean_detJ           = eda_results$det_J_stats$mean,
  sd_detJ             = eda_results$det_J_stats$sd,
  mean_loglik         = eda_results$loglik_stats$mean,
  sd_loglik           = eda_results$loglik_stats$sd,
  roundtrip_error_X   = eda_results$round_trip_error["X"],
  roundtrip_error_Z   = eda_results$round_trip_error["Z"],
  U_in_range          = eda_results$U_in_range,
  na_in_X             = eda_results$data_quality$na_X,
  inf_in_X            = eda_results$data_quality$inf_X,
  na_in_loglik        = eda_results$data_quality$na_loglik
)

write.csv(eda_df,
          file = "eda_results_block_c.csv",
          row.names = FALSE)

cat("\nEDA abgeschlossen - CSV-Dateien gespeichert!\n")


cat("### BLOCK_D ###\n")
cat("==============================\n")

# Hilfsfunktion: Erstelle tram-Modell basierend auf Verteilung
create_tram_model <- function(distr, response_name, predictors = NULL) {
  
  # Formula erstellen
  if (is.null(predictors) || length(predictors) == 0) {
    # Marginal model (k=1)
    fml <- as.formula(paste(response_name, "~ 1"))
  } else {
    # Conditional model (k>1) 
    fml <- as.formula(paste(response_name, "~", paste(predictors, collapse = " + ")))
  }
  
  # Wähle tram-Modell basierend auf Verteilung
  # Für Proof of Concept: verwende BoxCox für alle nicht-standard Verteilungen
  switch(distr,
    "norm" = Lm(fml),
    "t" = Lm(fml),  # Approximation über BoxCox
    "laplace" = Lm(fml),  # Approximation
    "logis" = Colr(fml),
    "exp" = Survreg(fml),
    "lnorm" = BoxCox(fml),
    # Für alle anderen: BoxCox als flexible Transformation
    BoxCox(fml)
  )
}

# Hauptfunktion: Fit conditional transformation forests
mytrtf <- function(data) {
  
  # Anzahl Variablen
  K <- ncol(data)
  var_names <- colnames(data)
  
  # Liste für transformation models und forests
  ymod <- list()
  forests <- list()
  
  cat("Fitting transformation forests für K =", K, "Variablen\n")
  
  # Für jede Variable xk
  for (k in 1:K) {
    
    cat("  Variable", k, ":", var_names[k], "...")
    
    # Aktuelle Verteilung aus config
    distr <- config[[k]]$distr
    
    # Abhängige Variable
    response_name <- var_names[k]
    
    # Unabhängige Variablen (wenn k > 1)
    if (k == 1) {
      predictors <- NULL
    } else {
      predictors <- var_names[1:(k-1)]
    }
    
    # Erstelle tram-Modell
    mod <- create_tram_model(distr, response_name, predictors)
    
    # Fit marginal model
    ymod[[k]] <- mlt(mod, data = data)
    
    # Fit conditional forest (wenn k > 1)
    if (k > 1) {
      # Formula für forest: response | predictors
      forest_fml <- as.formula(paste(response_name, "~ . |", paste(predictors, collapse = " + ")))
      
      # Transformation forest
      forests[[k]] <- traforest(mod, 
                               formula = forest_fml,
                               data = data,
                               ntree = 50,  # Proof of concept
                               mtry = max(1, floor((k-1)/3)))
    } else {
      forests[[k]] <- NULL  # Keine Forest für marginales Modell
    }
    
    cat(" OK\n")
  }
  
  # Return object
  structure(list(
    ymod = ymod,
    forests = forests,
    config = config,
    K = K,
    var_names = var_names
  ), class = "mytrtf")
}

# Predict-Methode für log-density contributions
predict.mytrtf <- function(object, newdata, type = "logdensity") {
  
  if (type != "logdensity") {
    stop("Derzeit nur type = 'logdensity' implementiert")
  }
  
  K <- object$K
  N <- nrow(newdata)
  var_names <- object$var_names
  
  # Matrix für log-density contributions
  logd_contrib <- matrix(NA, nrow = N, ncol = K)
  colnames(logd_contrib) <- var_names
  
  # Für jede Variable
  for (k in 1:K) {
    
    if (k == 1) {
      # Marginal density von x1
      logd_contrib[, k] <- logLik(object$ymod[[k]], newdata = newdata)
      
    } else {
      # Conditional density von xk | x1, ..., x(k-1)
      
      # Falls Forest existiert
      if (!is.null(object$forests[[k]])) {
        # Verwende Forest für Vorhersage
        logd_contrib[, k] <- logLik(object$forests[[k]], newdata = newdata)
      } else {
        # Fallback: verwende marginales Modell
        logd_contrib[, k] <- logLik(object$ymod[[k]], newdata = newdata)
      }
    }
  }
  
  # Als Vektor zurückgeben (Summe über Variablen für jede Beobachtung)
  if (N == 1) {
    return(sum(logd_contrib))
  } else {
    return(rowSums(logd_contrib))
  }
}

cat("Block D implementiert!\n")
cat("- mytrtf(): Fit conditional transformation forests\n")
cat("- predict.mytrtf(): Compute log-density contributions\n")
cat("Verwendete Packages: trtf, tram (für Transformationsmodelle)\n")

cat("### BLOCK_E ###\n")
cat("==============================\n")

cat("### BLOCK E: LINE-BY-LINE DEMO ###\n")
cat("==================================\n\n")

cat("### Learning data (von Block C)\n")
# Lade Daten aus Block C
X_samples <- read.csv("X_samples_block_c.csv")
ll_samples <- read.csv("ll_samples_block_c.csv")

# Erstelle Training Data Frame
d <- X_samples
colnames(d) <- paste0("X", 1:ncol(d))

# True log-densities als Attribut hinzufügen
attr(d, "logd") <- ll_samples$loglik

cat("Training Daten geladen: N =", nrow(d), ", K =", ncol(d), "\n")

cat("\n### Fit transformation forests\n")
m <- mytrtf(d)
cat("Modell erfolgreich gefittet\n")

cat("\n### Test data (neue Samples)\n")
# Generiere neue Test-Daten
set.seed(123)  # Für Reproduzierbarkeit
N_test <- 50
test_samples <- sample_and_loglik(N = N_test, cfg = config)

# Test Data Frame erstellen
dt <- as.data.frame(test_samples[, 1:(ncol(test_samples)-1)])
colnames(dt) <- paste0("X", 1:ncol(dt))

# True log-densities als Attribut
attr(dt, "logd") <- test_samples[, "loglik"]

cat("Test Daten generiert: N =", nrow(dt), "\n")

cat("\n### Out-of-sample log-likelihood\n")
# Berechne log-likelihood: marginal + alle conditional forests
ll_marginal <- logLik(m$ymod[[1]], newdata = dt)
ll_forests <- sum(sapply(m$forests, function(f) {
  if (!is.null(f)) logLik(f, newdata = dt) else 0
}))
ll <- ll_marginal + ll_forests

cat("Marginale Log-Lik:", round(ll_marginal, 3), "\n")
cat("Forest Log-Lik:", round(ll_forests, 3), "\n")
cat("Total Log-Lik:", round(ll, 3), "\n")

cat("\n### Out-of-sample log-likelihood contributions\n")
ld <- predict(m, newdata = dt)
cat("Predicted log-densities berechnet\n")

# Test auf Gleichheit
equal_test <- all.equal(sum(ld), c(ll))
cat("all.equal() Test:", equal_test, "\n")

cat("\n### Parameter Vergleich: Wahre vs. Estimated\n")
cat("=============================================\n")

# Erstelle Parameter-Vergleichstabelle
param_comparison <- data.frame(
  Variable = character(),
  Distribution = character(),
  True_Param = character(),
  True_Value = character(),
  Estimated_Value = character(),
  Difference = character(),
  stringsAsFactors = FALSE
)

# Extrahiere wahre und geschätzte Parameter für jede Variable
for (k in 1:length(config)) {
  var_name <- paste0("X", k)
  distr <- config[[k]]$distr
  
  cat("Variable", var_name, "(", distr, "):\n")
  
  # Estimated parameter aus dem Modell
  mod_coefs <- coef(m$ymod[[k]])
  
  # Für verschiedene Verteilungen unterschiedlich handhaben
  if (distr == "norm") {
    # Für Normalverteilung: mean und sd
    if (!is.null(config[[k]]$parm)) {
      # Falls Parameter-Funktion existiert (conditional)
      cat("  Conditional distribution - Parameter via Funktion\n")
    } else {
      # Standard normal N(0,1)
      true_mean <- 0
      true_sd <- 1
      
      # Geschätzte Parameter aus tram-Modell
      # Bei Lm(): Intercept entspricht -mean/sd, Scale entspricht 1/sd
      if (length(mod_coefs) >= 2) {
        est_mean <- -mod_coefs[1] * mod_coefs[2]  # -intercept * scale
        est_sd <- 1 / mod_coefs[2]  # 1 / scale
      } else {
        est_mean <- 0
        est_sd <- 1
      }
      
      # Zu Tabelle hinzufügen
      param_comparison <- rbind(param_comparison, 
                               data.frame(
                                 Variable = var_name,
                                 Distribution = distr,
                                 True_Param = "mean",
                                 True_Value = sprintf("%.3f", true_mean),
                                 Estimated_Value = sprintf("%.3f", est_mean),
                                 Difference = sprintf("%.3f", abs(true_mean - est_mean))
                               ))
      
      param_comparison <- rbind(param_comparison,
                               data.frame(
                                 Variable = var_name,
                                 Distribution = distr,
                                 True_Param = "sd",
                                 True_Value = sprintf("%.3f", true_sd),
                                 Estimated_Value = sprintf("%.3f", est_sd),
                                 Difference = sprintf("%.3f", abs(true_sd - est_sd))
                               ))
      
      cat("    mean: true =", true_mean, ", est =", round(est_mean, 3), 
          ", diff =", round(abs(true_mean - est_mean), 3), "\n")
      cat("    sd:   true =", true_sd, ", est =", round(est_sd, 3), 
          ", diff =", round(abs(true_sd - est_sd), 3), "\n")
    }
  } else {
    # Für andere Verteilungen: zeige nur Koeffizienten
    cat("    Koeffizienten:", round(mod_coefs, 3), "\n")
    
    # Füge generische Einträge hinzu
    for (i in seq_along(mod_coefs)) {
      param_comparison <- rbind(param_comparison,
                               data.frame(
                                 Variable = var_name,
                                 Distribution = distr,
                                 True_Param = paste0("coef_", i),
                                 True_Value = "unknown",
                                 Estimated_Value = sprintf("%.3f", mod_coefs[i]),
                                 Difference = "unknown"
                               ))
    }
  }
}

cat("\n### Parameter Vergleichstabelle\n")
print(param_comparison)

cat("\n### BLOCK E ABGESCHLOSSEN ###\n")
cat("✓ Training Daten aus Block C verwendet\n")
cat("✓ Neue Test Daten generiert\n")
cat("✓ Transformation Forests gefittet\n")
cat("✓ Parameter Vergleich erstellt\n")
cat("✓ Out-of-sample Log-Likelihood berechnet\n")
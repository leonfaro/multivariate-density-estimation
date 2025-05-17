suppressPackageStartupMessages({
  library(extraDistr)  # für Laplace
  library(tram)        # für Lm, Colr, BoxCox, mlt
  library(trtf)        # optional
  library(ggplot2)     # für Plots
})

# --- Block A: Konfiguration für K=4 ---
clip <- function(x, lo = 1e-6, hi = 1e6) pmin(hi, pmax(lo, x))
EPS <- 1e-10
config4 <- list(
  list(distr = "norm",    parm = NULL),
  list(distr = "t",       parm = function(d) list(df = 3 + 0.5 * d$X1)),
  list(distr = "laplace", parm = function(d) list(m = 0.3 * d$X2, s = clip(1 + 0.1 * d$X2))),
  list(distr = "logis",   parm = function(d) list(location = 0.2 * d$X3, scale = clip(1 + 0.05 * d$X3)))
)
config <- config4; K <- length(config)

# Numerisch stabile Hilfsfunktionen
dist_fun <- function(pref, name) get(paste0(pref, name))
safe_logpdf <- function(val, eps = EPS) log(pmax(val, eps))
safe_cdf    <- function(val, eps = EPS) pmax(eps, pmin(1 - eps, val))
get_pars <- function(k, x_prev, cfg) {
  ck <- cfg[[k]]
  if (is.null(ck$parm)) return(list())
  names(x_prev) <- paste0("X", seq_along(x_prev))
  ck$parm(as.data.frame(as.list(x_prev)))
}
pdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  pars <- get_pars(k, x_prev, cfg)
  dens <- do.call(dist_fun("d", cfg[[k]]$distr), c(list(x = xk), pars, list(log = FALSE)))
  if (log) safe_logpdf(dens) else dens
}
cdf_k <- function(k, xk, x_prev, cfg, log = TRUE) {
  dname <- cfg[[k]]$distr
  pars  <- get_pars(k, x_prev, cfg)
  args  <- c(
    if (dname == "sn") list(x = xk) else list(q = xk),
    pars,
    list(log.p = FALSE)
  )
  cdfv <- do.call(dist_fun("p", dname), args)
  cdfv <- safe_cdf(cdfv)
  if (log) log(cdfv) else cdfv
}
qtf_k <- function(k, u, x_prev, cfg, eps = EPS) {
  u <- pmin(1 - eps, pmax(eps, u)); dname <- cfg[[k]]$distr
  pars <- get_pars(k, x_prev, cfg)
  if (dname=="sn") {
    res <- tryCatch(
      do.call(dist_fun("q", dname), list(p = u, dp = c(unlist(pars), tau = 0), solver = "RFB")),
      error = function(e) do.call(dist_fun("q", dname), list(p = u, dp = c(unlist(pars), tau = 0), solver = "NR"))
    )
  } else {
    res <- do.call(dist_fun("q", dname), c(list(p = u), pars))
  }
  res[is.infinite(res)] <- sign(res[is.infinite(res)]) * 1e6
  res[is.na(res)]       <- 0
  res
}

# DGP via Triangular Transport
det_J <- function(logd) rowSums(logd)
loglik <- function(Z_eta, logdet) {
  -0.5 * rowSums(Z_eta^2) - (ncol(Z_eta)/2) * log(2*pi) + logdet
}
eta_sample <- function(N) {
  Z_eta <- matrix(rnorm(N * K), nrow = N)
  U_eta <- pnorm(Z_eta)
  list(U_eta = U_eta, Z_eta = Z_eta)
}
S_inv <- function(U_eta, cfg = config) {
  U_eta <- apply(U_eta, 2, function(u) pmin(1 - EPS, pmax(EPS, u)))
  n <- nrow(U_eta); X_pi <- matrix(NA_real_, n, K); logd <- matrix(NA_real_, n, K)
  for (i in seq_len(n)) {
    x_prev <- numeric(0)
    for (k in seq_len(K)) x_prev[k] <- qtf_k(k, U_eta[i,k], x_prev, cfg)
    X_pi[i,] <- x_prev
    for (k in seq_len(K)) {
      logd[i,k] <- pdf_k(k, X_pi[i,k], X_pi[i,seq_len(k-1)], cfg, log=TRUE) -
                   dnorm(qnorm(U_eta[i,k]), log=TRUE)
    }
  }
  list(X_pi = X_pi, U_eta = U_eta, Z_eta = qnorm(U_eta), logd = logd)
}
pi_sample <- function(N, cfg = config) {
  ref <- eta_sample(N); inv <- S_inv(ref$U_eta, cfg)
  list(X_pi = inv$X_pi, U_eta = ref$U_eta, Z_eta = inv$Z_eta, logd = inv$logd)
}

# --- Block C: EDA Diagnostics (Train, N_train >= 500) ---
set.seed(2044)
N_train <- 500
samp_train <- pi_sample(N_train)
X_pi_train <- samp_train$X_pi; U_eta_train <- samp_train$U_eta
Z_eta_train <- samp_train$Z_eta; logd_train <- samp_train$logd

detJ_train <- det_J(logd_train)
ll_train   <- loglik(Z_eta_train, detJ_train)
cat("Train EDA:\n")
cat("- X_pi_train Mean: ", paste(round(colMeans(X_pi_train),3), collapse=", "), "\n")
cat("- X_pi_train SD:   ", paste(round(apply(X_pi_train,2,sd),3), collapse=", "), "\n")
cat("- det(J)_train Range: [", round(min(detJ_train),3), ",", round(max(detJ_train),3), "]\n\n")

# --- Block C: EDA Diagnostics (Test, N_test >= 500) ---
set.seed(2045)
N_test <- 500
samp_test <- pi_sample(N_test)
X_pi_test <- samp_test$X_pi; U_eta_test <- samp_test$U_eta
Z_eta_test <- samp_test$Z_eta; logd_test <- samp_test$logd

detJ_test <- det_J(logd_test)
ll_test   <- loglik(Z_eta_test, detJ_test)
cat("Test EDA:\n")
cat("- X_pi_test Mean: ", paste(round(colMeans(X_pi_test),3), collapse=", "), "\n")
cat("- X_pi_test SD:   ", paste(round(apply(X_pi_test,2,sd),3), collapse=", "), "\n")
cat("- det(J)_test Range: [", round(min(detJ_test),3), ",", round(max(detJ_test),3), "]\n\n")

# Save CSVs für Debugging
if (!dir.exists("results")) dir.create("results")
write.csv(X_pi_train, "results/X_pi_train.csv", row.names=FALSE)
write.csv(U_eta_train, "results/U_eta_train.csv", row.names=FALSE)
write.csv(Z_eta_train, "results/Z_eta_train.csv", row.names=FALSE)
write.csv(logd_train,    "results/logd_train.csv",    row.names=FALSE)
write.csv(X_pi_test,  "results/X_pi_test.csv",  row.names=FALSE)
write.csv(U_eta_test, "results/U_eta_test.csv", row.names=FALSE)
write.csv(Z_eta_test, "results/Z_eta_test.csv", row.names=FALSE)
write.csv(logd_test,    "results/logd_test.csv",    row.names=FALSE)

# --- Block C: Parametrische Baseline (ohne Forest) ---
# MLE-Schätzung der Conditonal-Dists auf Trainingsdaten
nll_funs <- list(
  function(p,xs,Xprev) -sum(dnorm(xs, mean=p[1], sd=exp(p[2]), log=TRUE)),
  function(p,xs,Xprev) -sum(dt(xs, df=p[1]+p[2]*Xprev[,1], log=TRUE)),
  function(p,xs,Xprev) -sum(dlaplace(xs, m=p[1]*Xprev[,2], s=exp(p[2]*Xprev[,2]), log=TRUE)),
  function(p,xs,Xprev) -sum(dlogis(xs, location=p[1]*Xprev[,3], scale=exp(p[2]*Xprev[,3]), log=TRUE))
)
init_vals <- list(c(0,0), c(3,0.5), c(0.3,0.1), c(0.2,0.05))
param_est <- vector("list", K)
for(k in seq_len(K)) {
  xs    <- X_pi_train[,k]
  Xprev <- if(k>1) X_pi_train[,1:(k-1),drop=FALSE] else NULL
  fit   <- optim(init_vals[[k]], nll_funs[[k]], xs=xs, Xprev=Xprev)$par
  param_est[[k]] <- switch(k,
    list(mean=fit[1], sd=exp(fit[2])),
    list(a=fit[1], b=fit[2]),
    list(c=fit[1], e=fit[2]),
    list(g=fit[1], h=fit[2])
  )
}
# Out-of-Sample Loglikelihood auf Testdaten
true_ll_mat_test  <- sapply(seq_len(K), function(k)
  pdf_k(k, X_pi_test[,k], if(k>1) X_pi_test[,1:(k-1),drop=FALSE] else NULL, config, log=TRUE)
)
param_ll_mat_test <- sapply(seq_len(K), function(k) {
  xs    <- X_pi_test[,k]
  Xprev <- if(k>1) X_pi_test[,1:(k-1),drop=FALSE] else NULL
  p     <- param_est[[k]]
  switch(k,
    dnorm(xs, mean=p$mean, sd=p$sd, log=TRUE),
    dt(xs, df=p$a + p$b * Xprev[,1], log=TRUE),
    dlaplace(xs, m=p$c * Xprev[,2], s=exp(p$e * Xprev[,2]), log=TRUE),
    dlogis(xs, location=p$g * Xprev[,3], scale=exp(p$h * Xprev[,3]), log=TRUE)
  )
})
ll_delta_df_test <- data.frame(
  dim          = seq_len(K),
  distribution = sapply(config, `[[`, "distr"),
  ll_true_sum  = sapply(1:K, function(k) sum(true_ll_mat_test[,k])),
  ll_param_sum = sapply(1:K, function(k) sum(param_ll_mat_test[,k]))
)
ll_delta_df_test$delta_ll <- ll_delta_df_test$ll_true_sum - ll_delta_df_test$ll_param_sum
ll_delta_df_test[,3:5] <- round(ll_delta_df_test[,3:5],3)
print(ll_delta_df_test)

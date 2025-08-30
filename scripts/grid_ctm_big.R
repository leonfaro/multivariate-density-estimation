#!/usr/bin/env Rscript
# Larger grid search for Cross-term Map hyperparameters with N configurable

suppressMessages({
  source("main.R")
})

seed <- as.integer(Sys.getenv("SEED", 42))
set.seed(seed)
n <- as.integer(Sys.getenv("N", 100))
if (!is.finite(n) || n <= 0) n <- 100L
perm <- c(1L,2L,3L,4L)

# Bigger grid centered around prior best
Qs <- c(24L, 32L, 40L)
lambda_non_grid <- c(5e-3, 1e-2, 2e-2, 3e-2)
lambda_mon_grid <- c(3e-3, 5e-3, 7.5e-3, 1e-2, 2e-2)

dir.create("results", showWarnings = FALSE)
prep <- prepare_data(n, config, seed = seed)
S0 <- prep$S
S <- list(
  X_tr  = S0$X_tr[, perm, drop = FALSE],
  X_te  = S0$X_te[, perm, drop = FALSE]
)

# Baseline: Separable
sep <- trainSeparableMap(S, seed = seed)
sep_S <- sep$S
LD_sep <- -predict(sep_S, S$X_te, type = "logdensity_by_dim")
sum_sep <- sum(colMeans(LD_sep))
cat(sprintf("[BASELINE] Separable SUM NLL: %.3f with n=%d\n", sum_sep, n))

best <- list(nll = Inf, Q = NA_integer_, lambda_non = NA_real_, lambda_mon = NA_real_, se_sum = NA_real_)
rows <- list()

for (Q in Qs) for (l_non in lambda_non_grid) for (l_mon in lambda_mon_grid) {
  cat(sprintf("[BIG] Q=%d, lambda_non=%.4f, lambda_mon=%.4f ... ", Q, l_non, l_mon))
  t_start <- Sys.time()
  fit <- tryCatch(trainCrossTermMap(S, degree_g = 3, seed = seed,
                                    warmstart_from_separable = TRUE,
                                    sep_degree_g = 2, sep_lambda = 1e-3,
                                    Q = Q, lambda_non = l_non, lambda_mon = l_mon),
                  error = function(e) e)
  if (inherits(fit, "error")) { cat("ERROR\n"); 
    runtime_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    rows[[length(rows)+1]] <- data.frame(n=n, seed=seed, Q=Q, lambda_non=l_non, lambda_mon=l_mon,
                                         nll_sum_cross=NA_real_, se_sum_cross=NA_real_,
                                         nll_sum_sep=sum_sep, se_sum_sep=se_sum_sep,
                                         improved=FALSE, status="error", msg=as.character(fit$message),
                                         runtime=runtime_sec,
                                         stringsAsFactors = FALSE); next }
  LD_cross <- -predict(fit$S, S$X_te, type = "logdensity_by_dim")
  stopifnot(is.matrix(LD_cross), all(is.finite(LD_cross)))
  sum_cross <- sum(colMeans(LD_cross))
  improved <- sum_cross <= sum_sep
  runtime_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("SUM=%.3f %s\n", sum_cross, if (improved) "(<= separable)" else ""))
  rows[[length(rows)+1]] <- data.frame(n=n, seed=seed, Q=Q, lambda_non=l_non, lambda_mon=l_mon,
                                       nll_sum_cross=sum_cross,
                                       nll_sum_sep=sum_sep,
                                       improved=improved, status="ok", msg="",
                                       runtime=runtime_sec,
                                       stringsAsFactors = FALSE)
  if (sum_cross < best$nll) best <- list(nll=sum_cross, Q=Q, lambda_non=l_non, lambda_mon=l_mon)
}

df <- if (length(rows)>0) do.call(rbind, rows) else data.frame()
out_path <- sprintf("results/grid_ctm_big_seed%03d_n%03d.csv", seed, n)
if (nrow(df)>0) write.csv(df, out_path, row.names = FALSE)
cat(sprintf("[RESULTS] Saved big grid to %s\n", out_path))

cat("[BEST-BIG] Cross-term configuration by SUM NLL\n")
print(data.frame(Q = best$Q, lambda_non = best$lambda_non, lambda_mon = best$lambda_mon,
                 nll_sum_cross = best$nll,
                 nll_sum_sep = sum_sep,
                 improved = best$nll <= sum_sep,
                 stringsAsFactors = FALSE))

# Pretty ranking table: hyperparameters, nats (SUM NLL), runtime (s)
if (nrow(df) > 0) {
  tbl <- df
  # Keep only hyperparameters + metrics
  tbl <- tbl[, c("Q","lambda_non","lambda_mon","nll_sum_cross","runtime")]
  names(tbl) <- c("Q","lambda_non","lambda_mon","nats","runtime_sec")
  # Order by nats asc, then runtime
  ord <- order(tbl$nats, tbl$runtime_sec, na.last = TRUE)
  tbl <- tbl[ord, , drop = FALSE]
  rownames(tbl) <- NULL
  cat("[RANKED] Cross-term grid (best nats on top)\n")
  print(tbl)
}

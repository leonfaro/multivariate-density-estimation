source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ks_model.R")
source("models/ttm_model.R")
source("04_evaluation.R")

round_df <- function(df, digits = 3) {
  idx <- vapply(df, is.numeric, logical(1))
  df[idx] <- lapply(df[idx], round, digits = digits)
  df
}

N <- 50
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "beta", parm = function(d) list(shape1 = softplus(d$X2), shape2 = softplus(d$X1))),
  list(distr = "gamma", parm = function(d) list(shape = softplus(d$X3), scale = softplus(d$X2)))
)
perm <- c(3, 4, 1, 2)

#' @export
main <- function() {
  G <- list(
    N = N,
    config = config,
    seed = 42,
    split_ratio = 0.5
  )

  X <- gen_samples(G)                             # Script 2
  S <- train_test_split(X, G$split_ratio, G$seed) # Script 3

  ## Modelle in Originalreihenfolge
  M_TRUE <- fit_TRUE(S$X_tr, S$X_te, G$config)    # Script 5
  M_TRTF <- fit_TRTF(S$X_tr, S$X_te, G$config)    # Script models/trtf_model.R
  M_KS   <- fit_KS(S$X_tr, S$X_te, G$config)      # Script models/ks_model.R
  M_TTM  <- fit_TTM(S$X_tr, S$X_te)               # Script models/ttm_model.R

  baseline_ll <- logL_TRUE_dim(M_TRUE, S$X_te)
  trtf_ll     <- logL_TRTF_dim(M_TRTF, S$X_te)
  ks_ll       <- logL_KS_dim(M_KS, S$X_te)
  ttm_ll      <- logL_TTM_dim(M_TTM, S$X_te)

  tab_normal <- data.frame(
    dim = as.character(seq_along(G$config)),
    distribution = sapply(G$config, `[[`, "distr"),
    logL_baseline = baseline_ll,
    logL_trtf = trtf_ll,
    logL_ks = ks_ll,
    logL_ttm = ttm_ll,
    stringsAsFactors = FALSE
  )

  sum_row <- setNames(vector("list", ncol(tab_normal)), names(tab_normal))
  for (nm in names(tab_normal)) {
    if (nm == "dim") {
      sum_row[[nm]] <- "k"
    } else if (is.numeric(tab_normal[[nm]])) {
      sum_row[[nm]] <- sum(tab_normal[[nm]])
    } else {
      sum_row[[nm]] <- NA
    }
  }
  tab_normal <- rbind(tab_normal, as.data.frame(sum_row, stringsAsFactors = FALSE))

  ## Modelle mit Permutation
  X_tr_p <- S$X_tr[, perm, drop = FALSE]
  X_te_p <- S$X_te[, perm, drop = FALSE]
  colnames(X_tr_p) <- paste0("X", seq_len(ncol(X_tr_p)))
  colnames(X_te_p) <- paste0("X", seq_len(ncol(X_te_p)))

  M_TRUE_p <- fit_TRUE(X_tr_p, X_te_p, G$config)
  M_TRTF_p <- fit_TRTF(X_tr_p, X_te_p, G$config)
  M_KS_p   <- fit_KS(X_tr_p, X_te_p, G$config)
  M_TTM_p  <- fit_TTM(X_tr_p, X_te_p)

  baseline_ll_p <- logL_TRUE_dim(M_TRUE_p, X_te_p)
  trtf_ll_p     <- logL_TRTF_dim(M_TRTF_p, X_te_p)
  ks_ll_p       <- logL_KS_dim(M_KS_p, X_te_p)
  ttm_ll_p      <- logL_TTM_dim(M_TTM_p, X_te_p)

  tab_perm <- data.frame(
    dim = as.character(seq_along(G$config)),
    distribution = sapply(G$config, `[[`, "distr"),
    logL_baseline = baseline_ll_p,
    logL_trtf = trtf_ll_p,
    logL_ks = ks_ll_p,
    logL_ttm = ttm_ll_p,
    stringsAsFactors = FALSE
  )

  sum_row <- setNames(vector("list", ncol(tab_perm)), names(tab_perm))
  for (nm in names(tab_perm)) {
    if (nm == "dim") {
      sum_row[[nm]] <- "k"
    } else if (is.numeric(tab_perm[[nm]])) {
      sum_row[[nm]] <- sum(tab_perm[[nm]])
    } else {
      sum_row[[nm]] <- NA
    }
  }
  tab_perm <- rbind(tab_perm, as.data.frame(sum_row, stringsAsFactors = FALSE))

  ## Log-Dichten fuer Plot (Originalreihenfolge)
  ld_base <- rowSums(vapply(seq_along(G$config), function(k) {
    .log_density_vec(S$X_te[, k], G$config[[k]]$distr, M_TRUE$theta[[k]])
  }, numeric(nrow(S$X_te))))
  ld_trtf <- as.numeric(predict(M_TRTF, S$X_te, type = "logdensity"))
  ld_ks   <- as.numeric(predict(M_KS, S$X_te, type = "logdensity"))
  ld_ttm  <- as.numeric(predict(M_TTM, S$X_te, type = "logdensity"))

  plot(ld_base, ld_base, pch = 16, col = "black",
       xlab = "True log-Dichte",
       ylab = "Geschaetzte log-Dichte",
       main = "Originalreihenfolge")
  points(ld_base, ld_trtf, col = "blue", pch = 1)
  points(ld_base, ld_ks,   col = "red",  pch = 2)
  points(ld_base, ld_ttm,  col = "darkgreen", pch = 3)
  abline(a = 0, b = 1)
  legend("topleft", legend = c("true_baseline", "trtf", "ks", "ttm"),
         col = c("black", "blue", "red", "darkgreen"), pch = c(16, 1, 2, 3))

  ## Log-Dichten fuer Plot (Permutation)
  ld_base_p <- rowSums(vapply(seq_along(G$config), function(k) {
    .log_density_vec(X_te_p[, k], G$config[[k]]$distr, M_TRUE_p$theta[[k]])
  }, numeric(nrow(X_te_p))))
  ld_trtf_p <- as.numeric(predict(M_TRTF_p, X_te_p, type = "logdensity"))
  ld_ks_p   <- as.numeric(predict(M_KS_p,   X_te_p, type = "logdensity"))
  ld_ttm_p  <- as.numeric(predict(M_TTM_p,  X_te_p, type = "logdensity"))

  plot(ld_base_p, ld_base_p, pch = 16, col = "black",
       xlab = "True log-Dichte",
       ylab = "Geschaetzte log-Dichte",
       main = "Permutation")
  points(ld_base_p, ld_trtf_p, col = "blue", pch = 1)
  points(ld_base_p, ld_ks_p,   col = "red",  pch = 2)
  points(ld_base_p, ld_ttm_p,  col = "darkgreen", pch = 3)
  abline(a = 0, b = 1)
  legend("topleft", legend = c("true_baseline", "trtf", "ks", "ttm"),
         col = c("black", "blue", "red", "darkgreen"), pch = c(16, 1, 2, 3))

  print(round_df(tab_normal, digits = 3))
  print(round_df(tab_perm, digits = 3))

  invisible(list(normal = tab_normal, permutation = tab_perm))
}

# einfache Möglichkeit, weitere Modelle anzuhängen:
# - Quellcode in models/<neues>.R mit
#     fit_<NAME>()  und  logL_<NAME>()
# - anschließend hier laden und der Liste `models` hinzufügen.

if (sys.nframe() == 0L) {
  main()
}


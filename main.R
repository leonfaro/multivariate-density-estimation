source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ks_model.R")
source("models/ttm_model.R")
source("04_evaluation.R")
source("EDA.R")

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
  pdf("quick_test.pdf", width = 7, height = 7)
  on.exit(dev.off())
  t0 <- proc.time()

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

  t_true  <- system.time(
    baseline_ll <- logL_TRUE_dim(M_TRUE, S$X_te)
  )["elapsed"]
  t_trtf <- system.time(
    trtf_ll     <- logL_TRTF_dim(M_TRTF, S$X_te)
  )["elapsed"]
  t_ks   <- system.time(
    ks_ll       <- logL_KS_dim(M_KS, S$X_te)
  )["elapsed"]
  t_ttm  <- system.time(
    ttm_ll      <- logL_TTM_dim(M_TTM, S$X_te)
  )["elapsed"]

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

  runtime <- round((proc.time() - t0)[["elapsed"]], 2)
  message(paste0("Laufzeit: ", runtime, " Sekunden"))




  ## Log-Dichten fuer Plot (Permutation)
  ld_base_p <- rowSums(vapply(seq_along(G$config), function(k) {
    .log_density_vec(X_te_p[, k], G$config[[k]]$distr, M_TRUE_p$theta[[k]])
  }, numeric(nrow(X_te_p))))
  ld_trtf_p <- as.numeric(predict(M_TRTF_p, X_te_p, type = "logdensity"))
  ld_ks_p   <- as.numeric(predict(M_KS_p,   X_te_p, type = "logdensity"))
  ld_ttm_p  <- as.numeric(predict(M_TTM_p,  X_te_p, type = "logdensity"))

  scatter_data <- list(
    ld_base   = ld_base,
    ld_trtf   = ld_trtf,
    ld_ks     = ld_ks,
    ld_ttm    = ld_ttm,
    ld_base_p = ld_base_p,
    ld_trtf_p = ld_trtf_p,
    ld_ks_p   = ld_ks_p,
    ld_ttm_p  = ld_ttm_p
  )


  vec_normal <- c(true = t_true, trtf = t_trtf, ks = t_ks, ttm = t_ttm)
  kbl_tab <- combine_logL_tables(tab_normal, tab_perm,
                                 M_TRUE_p, M_TRTF_p, M_KS_p, M_TTM_p, X_te_p,
                                 vec_normal)

  create_EDA_report(S$X_tr, G$config,
                    scatter_data = scatter_data,
                    table_kbl = kbl_tab)


  invisible(kbl_tab)
}



if (sys.nframe() == 0L) {
  main()
}


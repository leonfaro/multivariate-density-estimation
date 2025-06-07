old_wd <- setwd("../..")
source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ks_model.R")
source("04_evaluation.R")

set.seed(1)
G <- setup_global()
G$N <- 10
X <- gen_samples(G)
S <- train_test_split(X, G$split_ratio, G$seed)

M_TRUE <- fit_TRUE(S$X_tr, S$X_te, G$config)
M_TRTF <- fit_TRTF(S$X_tr, S$X_te, G$config)
M_KS   <- fit_KS(S$X_tr, S$X_te, G$config)

t_true  <- system.time(logL_TRUE_dim(M_TRUE, S$X_te))["elapsed"]
t_trtf <- system.time(logL_TRTF_dim(M_TRTF, S$X_te))["elapsed"]
t_ks   <- system.time(logL_KS_dim(M_KS,   S$X_te))["elapsed"]

baseline_ll <- logL_TRUE_dim(M_TRUE, S$X_te)
trtf_ll     <- logL_TRTF_dim(M_TRTF, S$X_te)
ks_ll       <- logL_KS_dim(M_KS, S$X_te)

tab_normal <- data.frame(
  dim = as.character(seq_along(G$config)),
  distribution = sapply(G$config, `[[`, "distr"),
  logL_baseline = baseline_ll,
  logL_trtf = trtf_ll,
  logL_ks = ks_ll,
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

X_tr_p <- S$X_tr[, c(2,1,3,4)[seq_len(ncol(S$X_tr))], drop = FALSE]
X_te_p <- S$X_te[, c(2,1,3,4)[seq_len(ncol(S$X_te))], drop = FALSE]
colnames(X_tr_p) <- paste0("X", seq_len(ncol(X_tr_p)))
colnames(X_te_p) <- paste0("X", seq_len(ncol(X_te_p)))

M_TRUE_p <- fit_TRUE(X_tr_p, X_te_p, G$config)
M_TRTF_p <- fit_TRTF(X_tr_p, X_te_p, G$config)
M_KS_p   <- fit_KS(X_tr_p, X_te_p, G$config)

baseline_ll_p <- logL_TRUE_dim(M_TRUE_p, X_te_p)
trtf_ll_p     <- logL_TRTF_dim(M_TRTF_p, X_te_p)
ks_ll_p       <- logL_KS_dim(M_KS_p, X_te_p)

tab_perm <- data.frame(
  dim = as.character(seq_along(G$config)),
  distribution = sapply(G$config, `[[`, "distr"),
  logL_baseline = baseline_ll_p,
  logL_trtf = trtf_ll_p,
  logL_ks = ks_ll_p,
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

vec_normal <- c(true = t_true, trtf = t_trtf, ks = t_ks)

kbl_obj <- combine_logL_tables(tab_normal, tab_perm,
                              M_TRUE_p, M_TRTF_p, M_KS_p, NULL, X_te_p,
                              vec_normal)

setwd(old_wd)

test_that("combine_logL_tables returns kable", {
  expect_s3_class(kbl_obj, "knitr_kable")
})

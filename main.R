source("00_globals.R")
source("01_data_generation.R")
source("02_split.R")
source("models/true_model.R")
source("models/trtf_model.R")
source("models/ks_model.R")
source("models/ttm_base.R")
source("models/ttm_marginal.R")
source("04_evaluation.R")
source("replicate_code.R")

n <- 50
config <- list(
  list(distr = "norm", parm = NULL),
  list(distr = "exp",  parm = function(d) list(rate = softplus(d$X1))),
  list(distr = "beta",
       parm = function(d) list(shape1 = softplus(d$X2),
                               shape2 = softplus(d$X1))),
  list(distr = "gamma",
       parm = function(d) list(shape = softplus(d$X3),
                               scale = softplus(d$X2)))
)

#' Starte die komplette Analyse
#' 
#' Es werden zwei Log-Likelihood-Tabellen ausgegeben: einmal für die
#' normale Reihenfolge der Dimensionen und einmal für eine vorgegebene
#' Permutation.
#' 
#' @export
main <- function() {
  prep <- prepare_data(n, config, seed = 42)
  mods <- list(
    true = fit_TRUE(prep$S, config),
    trtf = fit_TRTF(prep$S, config),
    ks   = fit_KS(prep$S, config)
  )
  X_te <- prep$S$X_te
  res <- list()

  ll_true <- matrix(NA_real_, nrow = nrow(X_te), ncol = length(config))
  for (k in seq_along(config)) {
    ll_vec <- .log_density_vec(X_te[, k], config[[k]]$distr,
                               mods$true$theta[[k]])
    ll_true[, k] <- -ll_vec
  }
  res[["true"]] <- list(mean = colMeans(ll_true),
                         se = apply(ll_true, 2, sd)/sqrt(nrow(ll_true)))

  ll_trtf <- -predict(mods$trtf, X_te, type = "logdensity_by_dim")
  res[["trtf"]] <- list(mean = colMeans(ll_trtf),
                         se = apply(ll_trtf, 2, sd)/sqrt(nrow(ll_trtf)))

  ll_ks <- -predict(mods$ks, X_te, type = "logdensity_by_dim")
  res[["ks"]] <- list(mean = colMeans(ll_ks),
                      se = apply(ll_ks, 2, sd)/sqrt(nrow(ll_ks)))


  results_table <<- t(sapply(res, `[[`, "mean"))
  colnames(results_table) <<- paste0("dim", seq_len(length(config)))

  tab <- calc_loglik_tables(mods, config, X_te)
  print(tab)
  invisible(tab)
}

if (sys.nframe() == 0L) {
  main()
  replicate_code_scripts("main.R", "replicated_code.txt", env = globalenv())
}


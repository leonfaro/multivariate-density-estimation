#' Masked Autoregressive Flow – 5 Flows
fit_MAF <- function(S, config,
                    n_hidden = 100,
                    n_flows  = 5,
                    lr       = 1e-3,
                    n_epoch  = 200) {

  library(tensorflow); library(tfprobability)

  X_tr <- S$X_tr
  d    <- ncol(X_tr); n_batch <- 128

  base_dist <- tfd_multivariate_normal_diag(
    loc = rep(0, d), scale_diag = rep(1, d))

  maf_layer <- function()
    tfb_masked_autoregressive_flow(
      shift_and_log_scale_fn =
        tfb_masked_autoregressive_default_template(
          hidden_layers = list(n_hidden, n_hidden),
          activation = tf$nn$tanh))

  bijectors <- unlist(
    lapply(seq_len(n_flows),
           \(.) list(maf_layer(), tfb_permute(rev(seq_len(d)) - 1))),
    recursive = FALSE)

  flow_bij <- tfb_chain(rev(bijectors))
  dist     <- tfd_transformed_distribution(base_dist, flow_bij)

  opt <- tf$keras$optimizers$Adam(lr)
  ds  <- tensor_slices_dataset(X_tr) %>% dataset_shuffle(1000) %>%
    dataset_batch(n_batch)

  for (epoch in seq_len(n_epoch)) {
    itr <- as_iterator(ds)
    until_out_of_range({
      x <- iter_next(itr)
      with(tf$GradientTape() %as% tape, {
        loss <- -tf$reduce_mean(dist$log_prob(x))
      })
      grads <- tape$gradient(loss, dist$trainable_variables)
      opt$apply_gradients(
        purrr::transpose(list(grads, dist$trainable_variables)))
    })
  }

  structure(list(dist = dist), class = "maf_model")
}

predict.maf_model <- function(object, newdata,
                              type = c("logdensity",
                                       "logdensity_by_dim")) {
  type <- match.arg(type)
  lp   <- object$dist$log_prob(newdata) %>% as.numeric()
  if (type == "logdensity") return(lp)
  K <- ncol(newdata)
  matrix(rep(lp, each = K) / K, ncol = K, byrow = TRUE)
}

library(testthat)
source('../../scripts/halfmoon_plot.R')

dummy_model <- list(mu = c(0,0), sd = c(1,1))
class(dummy_model) <- 'dummy_model'
assign('predict.dummy_model', function(object, newdata, type = c('logdensity', 'logdensity_by_dim'), ...) {
  type <- match.arg(type)
  ld1 <- dnorm(newdata[,1], object$mu[1], object$sd[1], log=TRUE)
  ld2 <- dnorm(newdata[,2], object$mu[2], object$sd[2], log=TRUE)
  if (type == 'logdensity_by_dim') cbind(ld1, ld2) else ld1 + ld2
}, envir = .GlobalEnv)

test_that('eval_density_grid liefert korrekte Dimensionen und Cache', {
  xlim <- c(-1,1); ylim <- c(-1,1); gs <- 5L
  xseq <- seq(xlim[1], xlim[2], length.out = gs)
  yseq <- seq(ylim[1], ylim[2], length.out = gs)
  G <- as.matrix(expand.grid(xseq, yseq))
  res1 <- eval_density_grid(dummy_model, G, xlim, ylim, gs, seed=1, chunk=10, cores=1)
  expect_equal(dim(res1$LD), c(gs^2, 2))
  expect_equal(length(res1$joint), gs^2)
  expect_true(all(is.finite(res1$LD)))
  res2 <- eval_density_grid(dummy_model, G, xlim, ylim, gs, seed=1, chunk=10, cores=1)
  expect_identical(res1$joint, res2$joint)
})

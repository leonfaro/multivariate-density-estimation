
pkgs <- c("trtf", "kableExtra", "dplyr",
          "tibble", "tidyr")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}
library(parallel)
NC <- detectCores()
options(mc.cores = NC)

softplus <- function(x) log1p(exp(x))



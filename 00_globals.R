
pkgs <- c("trtf", "kableExtra", "dplyr",
          "tibble", "tidyr")

for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    library(p, character.only = TRUE)
  } else {
    message("Package ", p, " not installed; skipping.")
  }
}
if (requireNamespace("reticulate", quietly = TRUE))
  reticulate::use_python(Sys.which("python3"), required = TRUE)
library(parallel)
if (!exists("NC")) NC <- detectCores()
options(mc.cores = NC)

softplus <- function(x) log1p(exp(x))



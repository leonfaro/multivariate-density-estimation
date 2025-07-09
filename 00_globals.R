pkgs <- c("trtf", "dplyr", "tibble", "tidyr")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

if (!reticulate::py_module_available("tensorflow_probability")) {
  message("→ installiere TensorFlow 2.19 + TFP 0.23 (≈ kompatibel zu Py 3.12)…")
  tensorflow::install_tensorflow(
    version         = "2.19.0",                      # läuft mit Python 3.12
    extra_packages  = "tensorflow-probability==0.23.0",
    pip_ignore_installed = TRUE)
}


library(parallel)
if (!exists("NC")) NC <- detectCores()
options(mc.cores = NC)

softplus <- function(x) log1p(exp(x))



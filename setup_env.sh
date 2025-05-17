#!/usr/bin/env bash
set -e

REQUIRED_PACKAGES=(trtf tram ggplot2 extraDistr)

# Install R if missing
if ! command -v R >/dev/null 2>&1; then
  echo "R not found. Installing base R..." >&2
  sudo apt-get update
  sudo apt-get install -y --no-install-recommends gnupg ca-certificates software-properties-common
  sudo apt-key adv --no-tty --keyserver keyserver.ubuntu.com --recv-keys 'E298A3A825C0D65DFD57CBB651716619E084DAB9'
  sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
  sudo apt-get update
  sudo apt-get install -y r-base
fi

# Install required R packages
missing_pkgs=$(Rscript -e 'args <- commandArgs(TRUE); req <- args; miss <- req[!(req %in% rownames(installed.packages()))]; cat(paste(miss, collapse=" "))' "${REQUIRED_PACKAGES[@]}")
if [ -n "$missing_pkgs" ]; then
  echo "Installing R packages: $missing_pkgs" >&2
  Rscript -e "install.packages(strsplit('$missing_pkgs', ' ')[[1]], repos='https://cloud.r-project.org')"
fi

echo "R environment setup complete."

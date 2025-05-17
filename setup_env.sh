#!/usr/bin/env bash
set -e

REQUIRED_PACKAGES=(trtf tram ggplot2 extraDistr)

# Ensure R is available
if ! command -v R >/dev/null 2>&1; then
  echo "R not found. Please install R before running this script." >&2
  exit 1
fi

# Install required R packages
if [ -f renv.lock ]; then
  Rscript -e "if(!require('renv', quietly=TRUE)) install.packages('renv', repos=NULL, type='source'); renv::restore(prompt=FALSE)"
fi

echo "R environment setup complete."

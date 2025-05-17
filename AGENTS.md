# Instructions for Contributors

When working in this repository, ensure that a suitable R environment is available. The following packages are required:

- `trtf`
- `tram`
- `ggplot2`
- `extraDistr`

To set up the environment on Ubuntu with the most recent R release, run:

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends gnupg ca-certificates software-properties-common
sudo apt-key adv --no-tty --keyserver keyserver.ubuntu.com --recv-keys 'E298A3A825C0D65DFD57CBB651716619E084DAB9'
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt-get update
sudo apt-get install -y r-base
Rscript -e "install.packages(c('trtf', 'tram', 'ggplot2', 'extraDistr'), repos='https://cloud.r-project.org')"
```

These packages are needed for the entire repository.

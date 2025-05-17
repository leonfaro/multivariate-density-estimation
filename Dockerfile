FROM ghcr.io/sourcegraph/universal:latest

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        r-base r-base-dev && \
    rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('extraDistr','tram','trtf','ggplot2'), repos='https://cloud.r-project.org')"

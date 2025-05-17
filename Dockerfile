FROM scratch
COPY binaries/R-4.3.2-linux-x86_64.tar.gz /tmp/
RUN tar -xzf /tmp/R-4.3.2-linux-x86_64.tar.gz -C /opt && \
    ln -s /opt/R-4.3.2/bin/Rscript /usr/local/bin/Rscript
ENV RENV_PATHS_LIBRARY="/opt/R-4.3.2/site-library"
COPY renv/ renv/
COPY renv.lock renv.lock
RUN /usr/local/bin/Rscript -e "renv::restore(lockfile='renv.lock', repos=NULL)"
WORKDIR /workspace
COPY . .
ENTRYPOINT ["/usr/local/bin/Rscript"]
CMD ["scripts/run_all.R"]

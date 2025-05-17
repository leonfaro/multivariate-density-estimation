FROM rocker/r-base:latest
COPY renv.lock renv.lock
COPY renv/ renv/
RUN R -e "if(!require('renv')) install.packages('renv', repos=NULL, type='source'); renv::restore(prompt=FALSE)"
COPY . /app
WORKDIR /app
ENTRYPOINT ["/app/entrypoint.sh"]

FROM rocker/r-ver:4.3.3

# install2.r ist im Rocker-Image enthalten und kompiliert parallel
RUN install2.r --error \
      extraDistr \
      tram \
      trtf \
      ggplot2

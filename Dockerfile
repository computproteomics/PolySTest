FROM rocker/shiny
LABEL maintainer="Veit Schwaemmle <veits@bmb.sdu.dk>"
LABEL description="Docker image of PolySTest implementation on top of shiny-server. The number of to-be-installed R packages requires patience when building this image."

#RUN  sudo mount -o ro,remount /sys && sudo  mount -o rw,remount /sys
#RUN echo N | tee /sys/module/overlay/parameters/metacopy
#RUN rm -rf /var/cache/apt/* /var/lib/apt/lists/* /tmp/* /var/tmp/*
#RUN apt-get clean && apt-get update && apt-get install -y apt-utils

RUN apt-get update && apt-get install -y libssl-dev liblzma-dev libbz2-dev libicu-dev && apt-get clean 


RUN R -e "install.packages('BiocManager', repos='http://cran.us.r-project.org'); \
  update.packages(ask=F); \
  BiocManager::install(c('dplyr','plotly'),ask=F)"
RUN R -e "library(BiocManager); BiocManager::install(c('matrixStats','fdrtool','parallel','qvalue','circlize','DT','UpSetR','heatmaply','gplots','shinyBS','shinydashboard','limma'\
),ask=F)"


RUN rm -rf /srv/shiny-server
RUN mkdir /srv/shiny-server
COPY *R  /srv/shiny-server/
COPY *csv  /srv/shiny-server/
#COPY *pdf  /srv/shiny-server/

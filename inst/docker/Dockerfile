FROM rocker/shiny:4.4.1
LABEL maintainer="Veit Schwaemmle <veits@bmb.sdu.dk>"
LABEL description="Docker image of PolySTest implementation on top of shiny-server. The number of to-be-installed R packages requires patience when building this image."

#RUN  sudo mount -o ro,remount /sys && sudo  mount -o rw,remount /sys
#RUN echo N | tee /sys/module/overlay/parameters/metacopy
#RUN rm -rf /var/cache/apt/* /var/lib/apt/lists/* /tmp/* /var/tmp/*
#RUN apt-get clean && apt-get update && apt-get install -y apt-utils

RUN apt-get -o Acquire::AllowInsecureRepositories=true -o Acquire::AllowDowngradeToInsecureRepositories=true update && apt-get install -y libssl-dev liblzma-dev libbz2-dev libicu-dev && apt-get clean


RUN R -e "update.packages(ask=F); install.packages('BiocManager',ask=F); \
  BiocManager::install(c('dplyr','plotly','devtools'),ask=F)"

RUN rm -rf /srv/shiny-server/*
RUN mkdir /tmp/polystest
COPY .  /tmp/polystest/
COPY inst/shiny  /srv/shiny-server
RUN R -e "BiocManager::install(c('S4Vectors','XVector','S4Arrays'))"
RUN R -e "BiocManager::install(c('SparseArray','shinyBS','shinydashboard','shinyjs'))"
RUN cd /tmp/polystest && R -e "devtools::install()"
RUN chmod a+x /tmp/polystest/inst/cmdline/*
RUN cp /tmp/polystest/inst/cmdline/* /usr/bin/

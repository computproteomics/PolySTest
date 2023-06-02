#!/usr/bin/env Rscript
conda <- Sys.getenv("CONDA_PREFIX")
setwd(paste0(conda,"/share/polystest"))
shiny::runApp(port=3838)

#!/bin/bash 
R -e 'library(PolySTest); filePath <- system.file("extdata", "example_data.csv", package = "mypackage"); shiny::runApp(filePath, port=3838)'


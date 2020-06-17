#!/bin/bash 
pathname="\"$(readlink -f ./)\""
pathname=$pathname
R -e "shiny::runApp($pathname, port=3838)"


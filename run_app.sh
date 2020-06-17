#!/bin/bash 
pathname="\"$(readlink -f ./)\""
R -e "shiny::runApp($pathname, port=3838)"


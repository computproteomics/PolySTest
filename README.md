[![Anaconda-Server Badge](https://anaconda.org/bioconda/polystest/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

**Welcome to the repository of PolySTest**
developed at the [Protein Research Group](http://www.sdu.dk/en/Om_SDU/Institutter_centre/Bmb_biokemi_og_molekylaer_biologi/Forskning/Forskningsgrupper/Protein.aspx)  
Department of Biochemistry and Molecular Biology  
[University of Southern Denmark](http://www.sdu.dk)  

## Citation
When using PolySTest, please cite our paper:  
Veit Schwämmle, Christina E Hagensen, Adelina Rogowska-Wrzesinska, and Ole N. Jensen, doi: https://doi.org/10.1074/mcp.RA119.001777

PolySTest can be run as shiny app on our server or using the docker version of the tool to avoid installation issues and computational bottlenecks.

A command-line version is available to run PolySTest as workflow component (see below).

## Shiny app

### Web service

You can use our web server http://computproteomics.bmb.sdu.dk:

http://computproteomics.bmb.sdu.dk/Apps/PolySTest

Be aware that the tool does allow only one user to run the background R calculations at a time. Therefore the app might become temporarily irresponsive. However, multiple sessions are separated and your data won't be shared between sessions or overwritten. 

### Implementation on own computer

The easiest option is to use the docker image or as bioconda package.

#### Docker

`docker pull veitveit/polystest`

`docker run -t -i -p 3838:3838 veitveit/polystest`

and access the shiny app through http://localhost:3838

#### Bioconda

Install the package

`conda install -c bioconda polystest`

`run_app.sh`

and access the shiny app through http://localhost:3838


### Installation
Download the files from this repository

Install the following R libraries in R:
`library(BiocManager)
 BiocManager::install(c('matrixStats','fdrtool','parallel','qvalue','circlize','DT','UpSetR','heatmaply','gplots','shinyBS','shinydashboard','limma'),ask=F)`

You can run the shiny app from the server.R or ui.R files using [Rstudio](http://rstudio.com) or run the app on a [shiny-server](https://www.rstudio.com/products/shiny/shiny-server/)

Alternatively, just (in Linux) call the script `run_app.sh` and access the app through  http://localhost:3838

*Be aware that you need to have all files and all necessary R libraries described in the Installation*


### Build and use Docker image
A Dockerfile has been created on the basis of an OpenSuse distribution. Copy the repository to a folder and carry out the following command to build the images (takes a while)

`docker build -t veitveit/polystest .`

You can run the image by

`docker run -t -i -p 3838:3838 veitveit/polystest`

and access the server through http://localhost:3838

## Command-line version
You can run PolySTest via the file _runPolySTestCLI.R_. You need to specify a parameter file (yaml) to set parameters and the input file. See _polystest.yml_ for an example and extensive description of the parameters

### Installation
The command-line version comes with all version of PolySTest.


## Contact
For software issues and general questions, please submit an issue.

## License
GPL-2 or higher

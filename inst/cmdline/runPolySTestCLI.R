#!/usr/bin/env Rscript

# load libraries
library(yaml)
library(shiny)

## reading helper functions
# need to change to the source path and back
currPath <- getwd()
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(
  file.arg.name, "",
  initial.options[grep(file.arg.name, initial.options)]
)
script.basename <- dirname(script.name)
cat(paste0("Getting R functions from files located in ", script.basename), "\n")
setwd(script.basename)
source("HelperFuncs.R")
source("StatTestUnpaired.R")
source("StatTestPaired.R")
setwd(currPath)

######## reading parameter file
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("You need to specify the parameter file, see polystest.yml
        in the main folder as example and descriptions of the parameter values.
       This function needs to be carried out in the folder where the
       PolySTest files are located", call. = FALSE)
}
parfile <- args[1]

## read parameter file, show error otherwise
if(file.exists(parfile)) {
  cat(paste("Reading", parfile))
} else {
  stop(paste("Parameter file", parfile, "does not exist!"))
}
pars <- yaml.load_file(parfile)

validate_parameters(pars)

NumReps <- pars$numreps
NumCond <- pars$numcond
normalization <- pars$normalization
isPaired <- pars$paired
refCond <- pars$refcond
ColQuant <- pars$firstquantcol
delim <- pars$delim
if (delim == "tab") {
  delim <- "\t"
}
decimal <- pars$decimal
is_header <- pars$header
rep_grouped <- pars$rep_grouped
outfile <- pars$outfile
if (file.exists(outfile)) {
  warning(paste0(
    "******** WARNING: The file ", outfile,
    " will be overwritten!!"
  ))
}

filename <- pars$csvfile
if (!file.exists(filename)) {
  stop("csvfile does not exist", call. = FALSE)
}
threads <- pars$threads
# setting environment variable (if available)
Sys.setenv(SHINY_THREADS = threads)
Sys.getenv("SHINY_THREADS")

## Read file
dat <- read.csv(filename,
                header = is_header, sep = delim, dec = decimal,
                stringsAsFactors = F
)
rownames(dat) <- paste("feature", 1:nrow(dat))

# Data preparation
if (ColQuant - 1 > ncol(dat)) {
  stop("To large parameter firstquantcol!", call. = FALSE)
}
addInfo <- dat[, 1:(ColQuant - 1), drop = F]
for (c in 1:ncol(addInfo)) {
  addInfo[, c] <- as.character(addInfo[, c])
}
dat <- dat[, -(1:(ColQuant - 1))]

if (ncol(dat) < NumCond * NumReps) {
  stop("Not enough quantitative columns!", call. = F)
}

if (refCond > NumCond) {
  stop("Parameter refcond cannot be larger
        than the number of conditions", .call = F)
}

ndatcol <- ncol(dat)

# Reorder columns for replicates being grouped together
if (!rep_grouped) {
  print("reorder columns")
  act_cols <- rep(0:(NumCond - 1), NumReps) * NumReps +
    rep(1:(NumReps), each = NumCond)
  dat <- dat[, act_cols]
}
cat("\nExperimental setup as read:")
knitr::kable(matrix(colnames(dat),
                    ncol = NumCond, byrow = T,
                    dimnames = list(
                      rows = paste("Replicate", 1:NumReps),
                      cols = paste("Condition", 1:NumCond)
                    )
))

tncol <- 0
if (!is.null(addInfo)) {
  tncol <- ncol(dat) + ncol(addInfo)
} else {
  tncol <- ncol(dat)
}


# Normalize
if (normalization == "median") {
  dat <- t(t(dat) -
             colMedians(as.matrix(dat), na.rm = TRUE))
} else if (normalization == "mean") {
  dat <- t(t(dat) -
             colMeans(as.matrix(dat), na.rm = TRUE))
} else if (normalization == "sum_abs") {
  dat <- t(t(2^dat) /
             colSums(as.matrix(2^dat), na.rm = TRUE))
  dat <- log2(dat)
} else {
  dat <- limma::normalizeBetweenArrays(dat,
                                       method = normalization
  )
}

####### Defining conditions from names 
conditions <- update_conditions_with_lcs(dat, NumCond, NumReps, conditions = paste("C", 1:NumCond, sep = "")) 

####### Defining comparison for statistical tests
FullReg <- allComps <- NULL

if (NumCond > 1 & NumReps > 1) {
  allComps <- create_pairwise_comparisons(conditions, refCond)
  ncomps <- nrow(allComps)
  
  RR <- convert_comps_to_indices(allComps, conditions, NumCond, NumReps)

  # for maybe later to include onlyLIMMA option
  MAData <- NULL
  if (isPaired) {
    for (i in seq_len(ncol(RR))) {
      MAData <- cbind(MAData, dat[, RR[1, i]] - dat[, RR[2, i]])
    }
    rownames(MAData) <- rownames(dat)
    qvalues <- Paired(MAData, ncomps, NumReps)
  } else {
    qvalues <- UnpairedDesign(dat, RR, NumCond, NumReps)
  }
  MissingStats <- list(
    pNAvalues = matrix(1, ncol = ncomps, nrow = nrow(dat)),
    qNAvalues = matrix(1, ncol = ncomps, nrow = nrow(dat))
  )
  MissingStats <- MissingStatsDesign(dat, RR, NumCond, NumReps)
  
  # Prepare output data
  FullReg <- prepare_output_data(dat, qvalues, MissingStats, allComps)
  
} else {
  cat("!!!  Warning: Only one condition or one replicate per condition, no statistical tests performed  !!!\n")
  FullReg <- cbind(as.data.frame(addInfo), dat)
}



#### Save results
write.csv(FullReg, outfile)

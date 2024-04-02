#!/usr/bin/env Rscript

# load libraries
library(PolySTest)
library(SummarizedExperiment)
library(yaml)

## reading helper functions
# need to change to the source path and back
currPath <- getwd()
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(
  file.arg.name, "",
  initial.options[grep(file.arg.name, initial.options)]
)

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

PolySTest:::validate_parameters(pars)

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

# Data preparation
if (ColQuant - 1 > ncol(dat)) {
  stop("To large parameter firstquantcol!", call. = FALSE)
}

if (ncol(dat) < NumCond * NumReps) {
  stop("Not enough quantitative columns!", call. = F)
}

if (refCond > NumCond) {
  stop("Parameter refcond cannot be larger
        than the number of conditions", .call = F)
}

# Split data into feature metadata and quantitative data
addInfo <- dat[, seq_len(ColQuant - 1), drop = F]
dat <- dat[, -(seq_len(ColQuant - 1))]

# Reorder columns for replicates being grouped together
if (!rep_grouped) {
  print("reorder columns")
  act_cols <- rep(0:(NumCond - 1), NumReps) * NumReps +
    rep(1:(NumReps), each = NumCond)
  dat <- dat[, act_cols]
}

# Convert quantitative data to matrix if not already
quantDataMatrix <- as.matrix(dat)
rownames(quantDataMatrix) <- paste("feature", 1:nrow(dat))

# Create a DataFrame for sample metadata (colData)
# Assuming 'addInfo' contains feature-level metadata, adjust accordingly if it's sample-level metadata
sampleMetadata <- data.frame(Condition = rep(paste("Condition", 1:NumCond), NumReps),
                            Replicate = rep(1:NumReps, each=NumCond))

# Create the SummarizedExperiment object
fulldata <- SummarizedExperiment(assays = list(quant = quantDataMatrix),
                           colData = sampleMetadata)
rowData(fulldata) <- data.frame(addInfo)
# Adding metadata
metadata(fulldata) <- list(
  NumReps = NumReps,
  NumCond = NumCond
)

# Display experimental setup
cat("\nExperimental setup as read:")
print(knitr::kable(colData(fulldata)))

# Access the assay data
dat <- assay(fulldata, "quant")

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

# Update the assay data in 'fulldata'
assay(fulldata, "quant") <- dat

####### Defining conditions from names
fulldata <- update_conditions_with_lcs(fulldata)
conditions <- unique(colData(fulldata)$Condition)

####### Defining comparison for statistical tests
FullReg <- allComps <- NULL

if (NumCond > 1 & NumReps > 1) {
  allComps <- create_pairwise_comparisons(conditions, refCond)

  # Run tests via wrappers
  if (isPaired) {
    fulldata <- PolySTest_paired(fulldata, allComps)
  } else {
    fulldata <- PolySTest_unpaired(fulldata, allComps)
  }

} else {
  cat("!!!  Warning: Only one condition or one replicate per condition, no statistical tests performed  !!!\n")
}

FullReg <- cbind(rowData(fulldata), assay(fulldata))

#### Save results
write.csv(FullReg, outfile)

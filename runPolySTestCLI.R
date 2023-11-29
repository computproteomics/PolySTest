#!/usr/bin/env Rscript

library(yaml)
library(shiny)

######## reading parameter file
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("You need to specify the parameter file, see polystest.yml
        in the main folder as example and descriptions of the parameter values.
       This function needs to be carried out in the folder where the
       PolySTest files are located", call. = FALSE)
}
parfile <- args[1]
print(paste("Reading", parfile))

## read all parameters from yaml file + checks
pars <- yaml.load_file(parfile)
NumReps <- pars$numreps
if (!is.integer(NumReps)) {
  stop("Number of replicates is not a positive integer number", call. = FALSE)
}
if (NumReps < 0) {
  stop("Number of replicates is not a positive integer number", call. = FALSE)
}
NumCond <- pars$numcond
if (!is.integer(NumCond)) {
  stop("Number of conditions is not a positive integer number", call. = FALSE)
}
if (NumCond < 0) {
  stop("Number of conditions is not a positive integer number", call. = FALSE)
}
normalization <- pars$normalization
if (all(normalization != c(
  "none", "quantile", "median",
  "mean", "sum_abs", "cyclicloess"
))) {
  warning("No allowed normalization parameter value out of 'none',
  'quantile', 'median', 'mean', 'sum_abs' or 'cyclicloess'", call. = FALSE)
  normalization <- "none"
}
isPaired <- pars$paired
if (!is.logical(isPaired)) {
  stop("paired parameter should be 'true' or 'false'")
}
refCond <- pars$refcond
if (!is.integer(NumCond)) {
  stop("Parameter refCond is not a positive integer number", call. = FALSE)
}
if (NumCond < 0) {
  stop("Parameter refCond is not a positive integer number", call. = FALSE)
}
ColQuant <- pars$firstquantcol
if (!is.integer(ColQuant)) {
  stop("Parameter firstquantcol is not a positive integer number", call. = FALSE)
}
if (ColQuant < 1) {
  stop("Parameter firstquantcol is not a positive integer number > 0",
       call. = FALSE
  )
}
filename <- pars$csvfile
if (!file.exists(filename)) {
  stop("csvfile does not exist", call. = FALSE)
}
delim <- pars$delim
if (all(delim != c("tab", ";", ",", "|"))) {
  stop("No allowed delim parameter value out of ',', ';', '|' or 'tab'",
       call. = FALSE
  )
}
decimal <- pars$decimal
if (all(decimal != c(",", "."))) {
  stop("No allowed decimal parameter value out of ',' or '.'", call. = FALSE)
}
is_header <- pars$header
if (!is.logical(is_header)) {
  stop("header parameter should be 'true' or 'false'", call. = FALSE)
}
rep_grouped <- pars$rep_grouped
if (!is.logical(rep_grouped)) {
  stop("rep_grouped parameter should be 'true' or 'false'", call. = FALSE)
}
outfile <- pars$outfile
if (file.exists(outfile)) {
  warning(paste0(
    "******** WARNING: The file ", outfile,
    " will be overwritten!!"
  ))
}
threads <- pars$threads
if (!is.integer(threads)) {
  stop("Parameter threads is not a positive integer number", call. = FALSE)
}
if (threads < 1) {
  stop("Parameter threads is not a positive integer number > 0", call. = FALSE)
}
# setting environment variable
Sys.setenv(SHINY_THREADS = threads)
Sys.getenv("SHINY_THREADS")
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
setwd(currPath)

if (delim == "tab") {
  delim <- "\t"
}
dat <- read.csv(filename,
                header = is_header, sep = delim, dec = decimal,
                stringsAsFactors = F
)
rownames(dat) <- paste("feature", 1:nrow(dat))

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



conditions <- paste("C", 1:NumCond, sep = "")
## setting conditions names to common name of column names if there is any
# function to calculated "longest common substring" of 2 strings
# (from   https://stackoverflow.com/questions/35381180/identify-a-common-pattern)
lcs <- function(a, b) {
  A <- strsplit(a, "")[[1]]
  B <- strsplit(b, "")[[1]]
  L <- matrix(0, length(A), length(B))
  ones <- which(outer(A, B, "=="), arr.ind = TRUE)
  ones <- ones[order(ones[, 1]), , drop = F]
  if (nrow(ones) > 0) {
    for (i in seq_len(nrow(ones))) {
      v <- ones[i, , drop = FALSE]
      L[v] <- ifelse(any(v == 1), 1, L[v - 1] + 1)
    }
    paste0(A[(-max(L) + 1):0 + which(L == max(L),
                                     arr.ind = TRUE
    )[1]], collapse = "")
  } else {
    return(NA)
  }
}
for (i in 1:NumCond) {
  condnames <- colnames(dat)[0:(NumReps - 1) * NumCond + i]
  lcsname <- condnames[1]
  if (length(condnames) > 1) {
    for (j in 2:length(condnames)) {
      if (!is.na(lcsname)) {
        lcsname <- lcs(lcsname, condnames[j])
      }
    }
  }
  if (!is.na(lcsname) & sum(lcsname == conditions) == 0 ) {
    conditions[i] <- lcsname
  }
}

####### Defining comparison for statistical tests
# all vs. all or all vs. refCond?

FullReg <- allComps <- NULL
## Check for only one replicate per condition of only 1 condition
if (NumCond > 1 & NumReps > 1) {
  
  
  if (refCond > 0) {
    allComps <- cbind(conditions[refCond], conditions[conditions != conditions[refCond]])
  } else {
    for (i in seq_len(NumCond - 1)) {
      for (j in (i + 1):NumCond) {
        allComps <- rbind(allComps, c(conditions[i], conditions[j]))
      }
    }
  }
  
  cat("\nAll pairwise comparison between conditions:")
  colnames(allComps) <- paste("Condition", c("A", "B"))
  rownames(allComps) <- paste("Comparison", seq_len(nrow(allComps)))
  print(knitr::kable(t(allComps)))
  
  allComps <- unique(allComps)
  allComps <- allComps[allComps[, 1, drop = F] != allComps[, 2, drop = F], , drop = F]
  compNames <- paste(allComps[, 2], " vs ", allComps[, 1], sep = "")
  ncomps <- nrow(allComps)
  
  valComps <- matrix(NA, nrow = ncomps, ncol = ncol(allComps))
  for (i in seq_len(length(allComps))) {
    valComps[i] <-
      as.numeric(which(conditions == allComps[i]))
  }
  RR <- matrix(NA, ncol = ncomps * NumReps, nrow = 2)
  for (j in 1:ncomps) {
    compCond <- valComps[j, 2]
    refCond <- valComps[j, 1]
    RR[1, seq(j, ncomps * NumReps, ncomps)] <-
      seq(compCond, NumCond * NumReps, NumCond)
    RR[2, seq(j, ncomps * NumReps, ncomps)] <-
      seq(refCond, NumCond * NumReps, NumCond)
  }
  
  ###### Run tests
  
  # for maybe later to include onlyLIMMA option
  testNames <- c("limma", "Miss test", "rank products", "permutation test", "t-test")
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
  
  ### Preparing data for output
  LogRatios <- qvalues$lratios
  Pvalue <- cbind(
    qvalues$plvalues, MissingStats$pNAvalues, qvalues$pRPvalues,
    qvalues$pPermutvalues, qvalues$ptvalues
  )
  Qvalue <- cbind(
    qvalues$qlvalues, MissingStats$qNAvalues, qvalues$qRPvalues,
    qvalues$qPermutvalues, qvalues$qtvalues
  )
  Qvalue <- cbind(UnifyQvals(Qvalue, ncomps, 5), Qvalue)
  # print(head(LogRatios))
  
  colnames(LogRatios) <- paste("log-ratios", compNames)
  colnames(Pvalue) <- paste(
    "p-values", rep(testNames, each = ncomps),
    rep(compNames, length(testNames))
  )
  testNames2 <- c("PolySTest", testNames)
  colnames(Qvalue) <- paste(
    "FDR", rep(testNames2, each = ncomps),
    rep(compNames, length(testNames2))
  )
  fullData <- cbind(addInfo, dat)
  FullReg <- cbind(as.data.frame(fullData), LogRatios, Qvalue)
  
  #### report some statistics
  
  cat("------- summary of results --------\n")
  cat("Number of differentially regulated features with FDR < 0.01\n")
regMatr <- matrix(colSums(Qvalue < 0.01, na.rm = TRUE),
                    ncol = length(testNames2)
  )
  colnames(regMatr) <- testNames2
  rownames(regMatr) <- compNames
  print(knitr::kable(regMatr))
} else {
  cat("!!!  Warning: Only one condition or one replicate per condition, no statistical tests performed  !!!\n")
  FullReg <- cbind(as.data.frame(addInfo), dat)
}


#### Save results
write.csv(FullReg, outfile)

# Permutations: min 7 replicates, if not available, then generate from entire data

library(matrixStats)
library(fdrtool)
library(parallel)
library(qvalue)
library(limma)
source("rankprodbounds.R")

NumThreads <- 4
shiny_threads <- as.numeric(Sys.getenv("SHINY_THREADS"))
if (!is.na(shiny_threads)) {
  NumThreads <- shiny_threads
  print(paste("Set number of threads to", NumThreads))
}

# General parameters
NTests <- 1000 # for permutation tests
NumPermCols <- 7 # minimum number of columns for permutation tests (to ensure sufficient combinations)
NumRPPairs <- 100 # number of pairings in the unpaired rank product test to simulate an unpaired setting


#' Calculate statistics for permutation test
#'
#' This function calculates the statistics for a permutation test based on the input data.
#'
#' @param Data A matrix or data frame containing the data for the test.
#' @param Paired A logical value indicating whether the test is paired or not.
#' @return A numeric vector containing the calculated statistics.
#' @examples
#' Data <- matrix(rnorm(100), ncol = 10)
#' StatsForPermutTest(Data, Paired = FALSE)
#' @export
StatsForPermutTest <- function(Data, Paired) {
  if (Paired) {
    Stats <- rowMeans(Data, na.rm = TRUE) / rowSds(Data, na.rm = TRUE)
    Stats <- abs(Stats) * sqrt(rowSums(!is.na(Data)))
  } else {
    NumReps <- ncol(Data) * 0.5
    NumRealReps <- rowSums(!is.na(Data)) * 0.5
    Stats <- rowMeans(Data[, seq_len(NumReps)], na.rm = TRUE) -
      rowMeans(Data[, seq(NumReps + 1, length.out = NumReps)], na.rm = TRUE) /
        (sqrt(rowVars(Data[, seq_len(NumReps)], na.rm = TRUE) +
          rowVars(Data[, seq(NumReps + 1, length.out = NumReps)], na.rm = TRUE)))
    Stats <- abs(Stats) * NumRealReps
  }
  return(Stats)
}


#' Calculate the distribution of missing values
#'
#' This function calculates the distribution of missing values for a given number of repetitions and percentage of missing values.
#'
#' @param NumReps An integer indicating the number of repetitions.
#' @param PercNA A numeric value indicating the percentage of missing values.
#' @return A numeric vector containing the distribution of missing values.
#' @examples
#' MissValPDistr(10, 0.2)
#' @export
MissValPDistr <- function(NumReps, PercNA) {
  p <- PercNA
  d <- NumReps
  D <- rep(0, d + 1)
  # terms of binomial distribution
  binTerms <- NULL
  for (i in 0:d) {
    binTerms <- append(binTerms, choose(d, i) * p^i * (1 - p)^(d - i))
  }
  for (i in 0:d) {
    for (j in 0:(d - i)) {
      D[i + 1] <- D[i + 1] + binTerms[j + 1] * binTerms[j + i + 1]
    }
  }
  for (i in 1:d) {
    D[i + 1] <- 2 * D[i + 1]
  }
  D
}

#' RPStats Function
#'
#' This function calculates the p-values for the RP (Rank Product) statistic
#' based on the input tRPMAData and the number of replicates (NumReps).
#' The statistcs is one-sided, i.e. it only detects up-regulation.
#'
#' @param tRPMAData A matrix containing the expression data with rows as features and columns as replicates.
#' @param NumReps The number of replicates for each feature.
#'
#' @return A numeric vector containing the p-values for the RP statistic.
#'
#' @examples
#' tRPMAData <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
#' NumReps <- 3
#' RPStats(tRPMAData, NumReps)
#'
#' @export
RPStats <- function(tRPMAData, NumReps) {
  RPMADown_pvalues <- NULL
  NumElements <- rowSums(!is.na(tRPMAData))
  Rank <- NULL
  RP.own <- 0

  # Restrict to >0 replicated values
  iterNumEl <- unique(NumElements)
  iterNumEl <- iterNumEl[iterNumEl > 0]

  ## avoiding sets that are practically empty
  if (length(iterNumEl) == 0) {
    na_out <- as.numeric(rep(NA, nrow(tRPMAData)))
    names(na_out) <- rownames(tRPMAData)
    return(na_out)
  }

  # This function calculates p-values for RPMADown analysis based on the number of elements (d).
  # It returns a list of p-values for each iteration of d.
  RPMADown_pvalues <- lapply(iterNumEl, function(d) {
    tRPMADown_pvalues <- NULL
    # Subset the RPMAData based on the number of elements (d)
    RPMAData <- tRPMAData[NumElements == d, , drop = FALSE]
    # If d is greater than 1 and the number of columns in RPMAData is greater than the number of columns in tRPMAData
    if (d > 1 && length(as.matrix(RPMAData)) > ncol(tRPMAData)) {
      RP.own <- 0
      # Calculate the rank of each column in RPMAData and sum the ranks
      for (r in 1:NumReps) {
        Rank <- rank(RPMAData[, r], na.last = "keep") / (sum(!is.na(RPMAData[, r])) + 1)
        names(Rank) <- rownames(RPMAData)
        Rank[is.na(Rank)] <- 1
        RP.own <- RP.own + log(Rank)
      }
      # Calculate RP.own values and p-values
      RP.own <- exp(RP.own)
      tNumFeat <- length(RP.own)
      tRPout <- pgamma(-log(RP.own), d)
      names(tRPout) <- names(RP.own)
      tRPMADown_pvalues <- c(tRPMADown_pvalues, tRPout)
      # Return the calculated p-values
      tRPMADown_pvalues
    } else {
      # If d is less than or equal to 1 or the number of columns in RPMAData is not greater than the number of columns in tRPMAData,
      # return NA values for all rows in RPMAData
      na_out <- as.numeric(rep(NA, nrow(RPMAData)))
      names(na_out) <- rownames(RPMAData)
      return(na_out)
    }
  })
  return(unlist(RPMADown_pvalues))
}

#' Perform paired limma analysis
#'
#' This function performs paired limma analysis on MAData.
#'
#' @param MAData A matrix of gene expression data.
#' @param NumCond The number of conditions.
#' @param NumReps The number of replicates per condition.
#'
#' @return A list containing the p-values and q-values.
#' @keywords limma paired analysis
#' @export
limma_paired <- function(MAData, NumCond, NumReps) {
  MAReps <- rep(1:NumCond, NumReps)
  ## limma with ratios
  design <- plvalues <- NULL
  for (c in (1:(NumCond))) {
    design <- cbind(design, as.numeric(MAReps == c))
  }
  lm.fittedMA <- lmFit(MAData, design)
  lm.bayesMA <- eBayes(lm.fittedMA)
  topTable(lm.bayesMA)
  plvalues <- lm.bayesMA$p.value
  qlvalues <- matrix(NA, nrow = nrow(plvalues), ncol = ncol(plvalues), dimnames = dimnames(plvalues))
  # qvalue correction
  for (i in 1:ncol(plvalues)) {
    tqs <- qvalue(na.omit(plvalues[, i]))$qvalues
    qlvalues[names(tqs), i] <- tqs
  }
  return(plvalues = plvalues, qlvalues = qlvalues, Sds = sqrt(lm.bayesMA$s2.post))
}

# Perform paired t-tests
# This function performs paired t-tests on the given tMAData.
# Input: tMAData - a matrix of data for t-tests
# Output: a list containing the p-values (ptvalues) and q-values (qtvalues) of the t-tests
ttest_paired <- function(tMAData) {
  ## t-tests
  ptvalues <- sapply(1:nrow(tMAData), function(pep) {
    ifelse(sum(!is.na(tMAData[pep, ])) > 1,
      t.test(tMAData[pep, ])$p.value,
      NA
    )
  })
  names(ptvalues) <- rownames(tMAData)
  tqs <- qvalue(na.omit(ptvalues[, i]))$qvalues
  qtvalues <- rep(NA, length(ptvalues))
  names(qtvalues) <- names(ptvalues)
  qtvalues[names(tqs)] <- tqs

  return(list(ptvalues = ptvalues, qtvalues = qtvalues))
}

# Perform permutation tests
# This function performs permutation tests on the given tMAData.
# Input: tMAData - a matrix of data for permutation tests
# Output: a list containing the results of the permutation tests
permtest_paired <- function(tMAData) {
  ## Permutation tests
  # if necessary, add columns from randomized full set to reach min. NumPermCols replicates
  # randomizing also sign to avoid tendencies to one or the other side
  if (ncol(tMAData) < NumPermCols) {
    AddDat <- matrix(sample(as.vector(tMAData), (NumPermCols - ncol(tMAData)) * nrow(tMAData), replace = T), nrow = nrow(tMAData))
    PermMAData <- cbind(tMAData, AddDat)
  } else {
    PermMAData <- tMAData
  }

  RealStats <- StatsForPermutTest(tMAData, Paired = T)

  # run in parallel to speed up
  cl <- makeCluster(NumThreads)
  clusterExport(cl = cl, varlist = c("NumReps", "PermMAData", "RPStats", "StatsForPermutTest"), envir = environment())
  clusterEvalQ(cl = cl, library(matrixStats))

  PermutOut <- parLapply(cl, seq_len(NTests), function(x) {
    indat <- apply(
      PermMAData, 1,
      function(y) sample(y, NumReps) * sample(c(1, -1), NumReps, replace = T)
    )
    StatsForPermutTest(t(indat), T)
  })

  stopCluster(cl)

  PermutOut <- matrix(unlist(PermutOut), nrow = nrow(tMAData))
  pPermutvalues <- apply(cbind(RealStats, PermutOut), 1, function(x) ifelse(is.na(x[1]), NA, (1 + sum(x[1] < x[-1], na.rm = T)) / (sum(!is.na(x)))))
  qPermutvalues <- rep(NA, length(pPermutvalues))
  names(qPermutvalues) <- names(pPermutvalues)
  tqs <- p.adjust(na.omit(pPermutvalues[, i]), method = "BH")
  qPermutvalues[names(tqs), i] <- tqs

  return(list(pPermutvalues = pPermutvalues, qPermutvalues = qPermutvalues))
}

#' Perform unpaired limma analysis
#' This function performs unpaired limma analysis on Data.
#' Input: Data - a matrix of gene expression data
#'       NumCond - the number of conditions
#'      NumReps - the number of replicates per condition
#' Output: a list containing the p-values (plvalues) and q-values (qlvalues) of the limma analysis
#' @keywords limma unpaired analysis
#' @export
#' @examples
#' Data <- matrix(rnorm(100), ncol = 5)
#' results <- Unpaired(Data, 5, 2)
#' print(results)
#' @export
#' @keywords limma unpaired analysis
limma_unpaired <- function(Data, NumCond, NumReps, RRCateg) {
  Reps <- rep(1:NumCond, NumReps)
  NumComps <- ncol(RRCateg)
  design <- model.matrix(~ 0 + factor(Reps - 1))
  colnames(design) <- paste("i", c(1:NumCond), sep = "")
  contrasts <- NULL
  First <- 1
  for (i in (1:NumComps)) {
    contrasts <- append(contrasts, paste(colnames(design)[RRCateg[2, i]], "-", colnames(design)[RRCateg[1, i]], sep = ""))
  }
  contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)
  lm.fitted <- lmFit(Data, design)

  lm.contr <- contrasts.fit(lm.fitted, contrast.matrix)
  lm.bayes <- eBayes(lm.contr)
  topTable(lm.bayes)
  plvalues <- lm.bayes$p.value
  qlvalues <- matrix(NA, nrow = nrow(plvalues), ncol = ncol(plvalues), dimnames = dimnames(plvalues))
  # qvalue correction
  for (i in 1:ncol(plvalues)) {
    tqs <- qvalue(na.omit(plvalues[, i]))$qvalues
    qlvalues[names(tqs), i] <- tqs
  }

  return(list(plvalues = plvalues, qlvalues = qlvalues, Sds = sqrt(lm.bayes$s2.post)))
}

#' Perform unpaired t-tests on two datasets
#'
#' This function performs unpaired t-tests between corresponding rows of two datasets.
#' It calculates the p-values and q-values for each row, indicating the significance of the difference between the two datasets.
#'
#' @param tData The first dataset, a matrix or data frame
#' @param trefData The second dataset, a matrix or data frame
#'
#' @return A list containing the p-values and q-values for each row
#'
#' @examples
#' tData <- matrix(rnorm(100), nrow = 10)
#' trefData <- matrix(rnorm(100), nrow = 10)
#' result <- ttest_unpaired(tData, trefData)
#' print(result$ptvalues)
#' print(result$qtvalues)
#'
#' @export
ttest_unpaired <- function(tData, trefData) {
  ## t-tests
  tptvalues <- sapply(1:nrow(tData), function(pep) {
    ifelse(sum(!is.na(tData[pep, ])) > 1 & sum(!is.na(trefData[pep, ])) > 1,
      t.test(unlist(tData[pep, ]), unlist(trefData[pep, ]))$p.value,
      NA
    )
  })
  names(tptvalues) <- rownames(tData)
  ptvalues <- tptvalues
  tqs <- qvalue(na.omit(ptvalues))$qvalues
  qtvalues <- rep(NA, length(ptvalues))
  names(qtvalues) <- names(ptvalues)
  qtvalues[names(tqs)] <- tqs

  return(list(ptvalues = ptvalues, qtvalues = qtvalues))
}


#' rp_unpaired function
#'
#' This function calculates the p-values and q-values for unpaired random pairing combinations.
#'
#' @param tData The data matrix for the test group.
#' @param trefData The data matrix for the reference group.
#' @param NumReps The number of replicates.
#'
#' @return A list containing the p-values and q-values.
#'
#' @examples
#' tData <- matrix(rnorm(1000), nrow = 100)
#' trefData <- matrix(rnorm(1000), nrow = 100)
#' NumReps <- 10
#' rp_unpaired(tData, trefData, NumReps)
#'
#' @export
rp_unpaired <- function(tData, trefData, NumReps) {
  # calculate NumRPPairs random pairing combinations and then take mean of p-values
  tpRPvalues <- matrix(NA,
    ncol = NumRPPairs, nrow = nrow(tData),
    dimnames = list(rows = rownames(tData), cols = 1:NumRPPairs)
  )
  cl <- makeCluster(NumThreads)
  clusterExport(cl = cl, varlist = c("NumReps", "tData", "trefData", "RPStats"), envir = environment())
  clusterEvalQ(cl = cl, library(matrixStats))

  RPparOut <- parallel::parLapply(cl, 1:NumRPPairs, function(x) {
    tRPMAData <- tData[, sample(1:NumReps)] - trefData[, sample(1:NumReps)]
    # Up
    RPMAUp_pvalues <- RPStats(tRPMAData, NumReps)
    # Down
    RPMADown_pvalues <- RPStats(-tRPMAData, NumReps)
    ttt <- rowMins(cbind(RPMAUp_pvalues, RPMADown_pvalues), na.rm = T) * 2
    ttt[ttt > 1] <- 1
    names(ttt) <- names(RPMAUp_pvalues)
    ttt
  })
  stopCluster(cl)

  for (p in 1:NumRPPairs) {
    # print(RPparOut[[p]])
    tpRPvalues[names(RPparOut[[p]]), p] <- RPparOut[[p]]
  }
  tpRPvalues[!is.finite(tpRPvalues)] <- NA
  pRPvalues <- rowMeans(tpRPvalues, na.rm = T)
  qRPvalues <- rep(NA, length(pRPvalues))
  names(qRPvalues) <- names(pRPvalues)
  tqs <- p.adjust(na.omit(pRPvalues), method = "BH")
  qRPvalues[names(tqs)] <- tqs

  return(list(pRPvalues = pRPvalues, qRPvalues = qRPvalues))
}

#' perm_unpaired function
#'
#' This function performs permutation testing for unpaired data.
#'
#' @param tData The test data matrix.
#' @param trefData The reference data matrix.
#' @param NumReps The number of repetitions for permutation testing.
#'
#' @return A list containing the p-values and q-values for the permutation test.
#'
#' @details The function adds columns from the randomized full set to reach the minimum number of permutation columns (NumPermCols) replicates. It randomizes the sign as well to avoid tendencies to one or the other side. In the unpaired case, it also normalizes by the mean of the entire sample to avoid strange effects. The function then performs permutation testing using parallel computing, and calculates the p-values and q-values based on the permutation results.
#'
#' @examples
#' tData <- matrix(rnorm(1000), nrow = 100)
#' trefData <- matrix(rnorm(1000), nrow = 100)
#' NumReps <- 10
#' result <- perm_unpaired(tData, trefData, NumReps)
#'
#' @export
perm_unpaired <- function(tData, trefData, NumReps) {
  # add columns from randomized full set to reach min. NumPermCols replicates
  # randomizing also sign to avoid tendencies to one or the other side
  # In the unpaired case, also normalize by mean of the entire sample to avoid strange effects
  tData <- tData - mean(as.numeric(unlist(tData)), na.rm = T)
  trefData <- trefData - mean(as.numeric(unlist(trefData)), na.rm = T)
  if (ncol(tData) * 2 < NumPermCols) {
    AddDat <- matrix(sample(as.vector(unlist(tData)), (NumPermCols * 0.5 - ncol(tData)) * nrow(tData), replace = T), nrow = nrow(tData))
    PermData <- cbind(tData, AddDat)
    AddDat <- matrix(sample(as.vector(unlist(trefData)), (NumPermCols * 0.5 - ncol(trefData)) * nrow(trefData), replace = T), nrow = nrow(trefData))
    PermFullData <- cbind(PermData, trefData, AddDat)
  } else {
    PermFullData <- cbind(tData, trefData)
  }
  RealStats <- StatsForPermutTest(as.matrix(cbind(trefData, tData)), Paired = F)
  # print(head(PermFullData))

  cl <- makeCluster(NumThreads)
  clusterExport(cl = cl, varlist = c("NumReps", "PermFullData", "RPStats", "StatsForPermutTest"), envir = environment())
  clusterEvalQ(cl = cl, library(matrixStats))
  PermutOut <- parallel::parLapply(cl, 1:NTests, function(x) {
    indat <- apply(PermFullData, 1, function(y) sample(y, NumReps * 2) * sample(c(1, -1), NumReps * 2, replace = T))
    StatsForPermutTest(t(indat), F)
  })
  stopCluster(cl)

  PermutOut <- matrix(unlist(PermutOut), nrow = nrow(tData))
  PermutOut[!is.finite(PermutOut)] <- NA
  RealStats[!is.finite(RealStats)] <- NA
  pPermutvalues <- apply(cbind(RealStats, PermutOut), 1, function(x) ifelse(is.na(x[1]) | sum(!is.na(x)) == 0, NA, (1 + sum(x[1] < x[-1], na.rm = T)) / (sum(!is.na(x)))))
  qPermutvalues <- rep(NA, length(pPermutvalues))
  names(qPermutvalues) <- names(pPermutvalues)
  tqs <- qvalue(na.omit(pPermutvalues))$qvalues
  qPermutvalues[names(tqs)] <- tqs

  return(list(pPermutvalues = pPermutvalues, qPermutvalues = qPermutvalues))
}

#' Paired Test Analysis
#'
#' This function performs paired test analysis on the given MAData.
#'
#' @param MAData A matrix containing the data for the analysis.
#' @param NumCond The number of conditions in the analysis.
#' @param NumReps The number of replicates per condition.
#'
#' @return A list containing the results of the paired test analysis, including p-values, q-values, and ratios.
#'
#' @examples
#' MAData <- matrix(rnorm(100), ncol = 5)
#' results <- Paired(MAData, 5, 2)
#' print(results)
#'
#' @export
Paired <- function(MAData, NumCond, NumReps) {
  MAReps <- rep(1:NumCond, NumReps)

  cat("Running paired tests\n")
  ## limma with ratios
  limma_out <- limma_paired(MAData, NumCond, NumReps)
  plvalues <- limma_out$plvalues
  qlvalues <- limma_out$qlvalues
  Sds <- limma_out$Sds

  cat("limma completed\n")

  ## p-value objects
  ptvalues <- NULL
  pRPvalues <- matrix(NA,
    ncol = NumCond, nrow = nrow(MAData),
    dimnames = list(rows = rownames(MAData), cols = paste("RP p-values", 1:NumCond))
  )
  pPermutvalues <- matrix(NA,
    ncol = NumCond, nrow = nrow(MAData),
    dimnames = list(rows = rownames(MAData), cols = paste("Permutation p-values", 1:NumCond))
  )

  ## q-value objects
  qRPvalues <- qtvalues <- qPermutvalues <- matrix(NA, nrow = nrow(MAData), ncol = NumCond, dimnames = list(rows = rownames(MAData), cols = 1:NumCond))

  lratios <- NULL

  cat("Running rank products and permutations tests ...\n")
  pb <- txtProgressBar(0.9, NumCond)
  for (vs in 1:NumCond) {
    if (!is.null(getDefaultReactiveDomain())) {
      setProgress(0.1 + 0.3 / NumCond * vs, detail = paste("tests for comparison", vs, "of", NumCond))
    }

    tMAData <- MAData[, MAReps == vs]

    ## t-tests
    ttest_out <- ttest_paired(tMAData)
    ptvalues <- cbind(ptvalues, ttest_out$ptvalues)
    qtvalues[, vs] <- ttest_out$qtvalues

    ## rank products
    tRPMAData <- MAData[, MAReps == vs]
    # Up
    RPMAUp_pvalues <- RPStats(tRPMAData, NumReps)
    # Down
    RPMADown_pvalues <- RPStats(-tRPMAData, NumReps)
    ttt <- rowMins(cbind(RPMAUp_pvalues, RPMADown_pvalues), na.rm = T) * 2
    ttt[ttt > 1] <- 1
    pRPvalues[names(RPMAUp_pvalues), vs] <- ttt
    tqs <- p.adjust(na.omit(pRPvalues[, vs]), method = "BH")
    qRPvalues[names(tqs), vs] <- tqs

    ## Permutation tests
    perm_out <- permtest_paired(tMAData, NumReps)
    pPermutvalues[, vs] <- perm_out$pPermutvalues
    qPermutvalues[, vs] <- perm_out$qPermutvalues

    lratios <- cbind(lratios, rowMeans(MAData[, MAReps == i], na.rm = T))
    setTxtProgressBar(pb, vs)
  }
  cat("rank products and permutation test completed\n")
  close(pb)

  return(list(
    lratios = lratios, ptvalues = ptvalues, plvalues = plvalues,
    pRPvalues = pRPvalues, pPermutvalues = pPermutvalues,
    qtvalues = qtvalues, qlvalues = qlvalues,
    qRPvalues = qRPvalues, qPermutvalues = qPermutvalues,
    Sds = Sds))
}


#' Unpaired
#'
#' This function performs significance analysis for unpaired data.
#'
#' @param Data A matrix of gene expression data.
#' @param NumCond The number of conditions.
#' @param NumReps The number of replicates per condition.
#'
#' @return A list containing various results of the significance analysis, including lratios, ptvalues, plvalues, pRPvalues, pPermutvalues, qtvalues, qlvalues, qRPvalues, qPermutvalues, and Sds.
#'
#' @examples
#' Data <- matrix(rnorm(100), ncol = 5)
#' result <- Unpaired(Data, 5, 2)
#' print(result)
#'
#' @export
# for comparison versus first condition in table
Unpaired <- function(Data, NumCond, NumReps) {
  ##########################################################
  # significance analysis
  Reps <- rep(1:NumCond, NumReps)

  # Normalize row-wise by mean
  Data <- Data - rowMeans(Data, na.rm = T)


  ## limma
  cat("Running limma\n")
  limma_out <- limma_unpaired(Data, NumCond, NumReps, rbind(1, seq(2, NumCond, by = 1)))
  plvalues <- limma_out$plvalues
  qlvalues <- limma_out$qlvalues

  cat("limma completed\n")

  ## p-value objects
  ptvalues <- qtvalues <- matrix(NA, nrow = nrow(Data), ncol = NumCond - 1, dimnames = list(rows = rownames(Data), cols = 1:(NumCond - 1)))
  pRPvalues <- matrix(NA, ncol = NumCond - 1, nrow = nrow(Data), dimnames = list(rows = rownames(Data), cols = paste("RP p-values", 1:(NumCond - 1))))
  pPermutvalues <- matrix(NA, ncol = NumCond - 1, nrow = nrow(Data), dimnames = list(rows = rownames(Data), cols = paste("Permutation p-values", 1:(NumCond - 1))))
  ## q-value objects
  qRPvalues <- qtvalues <- qPermutvalues <- matrix(NA, nrow = nrow(Data), ncol = NumCond - 1, dimnames = list(rows = rownames(Data), cols = 1:(NumCond - 1)))

  lratios <- NULL

  cat("Running rank products and permutations tests ...\n")
  pb <- txtProgressBar(1.9, NumCond)
  for (vs in seq(from = 2, to = NumCond, by = 1)) {
    if (!is.null(getDefaultReactiveDomain())) {
      setProgress(0.1 + 0.3 / (NumCond - 1) * vs, detail = paste("tests for comparison", vs - 1, "of", NumCond - 1))
    }
    tData <- Data[, Reps == vs]
    trefData <- Data[, Reps == 1]


    ## t-test
    ttest_out <- ttest_unpaired(tData, trefData)
    ptvalues[, vs - 1] <- ttest_out$ptvalues
    qtvalues[, vs - 1] <- ttest_out$qtvalues

    ## rank products
    rp_out <- rp_unpaired(tData, trefData, NumReps)
    pRPvalues[, vs - 1] <- rp_out$pRPvalues
    qRPvalues[, vs - 1] <- rp_out$qRPvalues


    ## Permutation tests
    perm_out <- perm_unpaired(tData, trefData, NumReps)
    pPermutvalues[, vs - 1] <- perm_out$pPermutvalues
    qPermutvalues[, vs - 1] <- perm_out$qPermutvalues

    lratios <- cbind(lratios, rowMeans(Data[, Reps == vs], na.rm = T) - rowMeans(Data[, Reps == 1], na.rm = T))
    setTxtProgressBar(pb, vs)
  }
  close(pb)

  return(list(
    lratios = lratios, ptvalues = ptvalues, plvalues = plvalues, pRPvalues = pRPvalues, pPermutvalues = pPermutvalues,
    qtvalues = qtvalues, qlvalues = qlvalues, qRPvalues = qRPvalues, qPermutvalues = qPermutvalues, Sds = sqrt(lm.bayes$s2.post)
  ))
}

#' UnpairedDesign function
#'
#' This function performs significance analysis using different statistical tests such as limma, rank products, and permutation tests.
#'
#' @param Data A matrix of gene expression data.
#' @param RR A matrix specifying the experimental design.
#' @param NumCond The number of conditions in the experiment.
#' @param NumReps The number of replicates per condition.
#'
#' @return A list containing the following results:
#'   - lratios: The log ratios of gene expression between conditions.
#'   - ptvalues: The p-values from t-tests.
#'   - plvalues: The p-values from limma tests.
#'   - pRPvalues: The p-values from rank products tests.
#'   - pPermutvalues: The p-values from permutation tests.
#'   - qtvalues: The q-values from t-tests.
#'   - qlvalues: The q-values from limma tests.
#'   - qRPvalues: The q-values from rank products tests.
#'   - qPermutvalues: The q-values from permutation tests.
#'   - Sds: The standard deviations of the Bayesian linear model.
#'
#' @examples
#' # Example usage of UnpairedDesign function
#' data <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' RR <- matrix(c(1, 2, 2, 3), nrow = 2, ncol = 2)
#' result <- UnpairedDesign(data, RR, 3, 5)
#' print(result)
#'
#' @export
UnpairedDesign <- function(Data, RR, NumCond, NumReps) {
  ##########################################################
  # significance analysis
  Reps <- rep(1:NumCond, NumReps)

  # Normalize row-wise by mean
  Data <- Data - rowMeans(Data, na.rm = T)

  # Number of tests
  NumComps <- ncol(RR) / NumReps
  RRCateg <- RR[, 1:NumComps, drop = F]

  ## limma
  cat("Running limma tests\n")
  lm_out <- limma_unpaired(Data, NumCond, NumReps, RRCateg)
  plvalues <- lm_out$plvalues
  qlvalues <- lm_out$qlvalues
  Sds <- lm_out$Sds

  ## rank products + t-test
  ptvalues <- qtvalues <- matrix(NA, nrow = nrow(Data), ncol = NumComps, dimnames = list(rows = rownames(Data), cols = 1:NumComps))
  pRPvalues <- qRPvalues <- matrix(NA, ncol = NumComps, nrow = nrow(Data), dimnames = list(rows = rownames(Data), cols = paste("RP p-values", 1:(NumComps))))
  pPermutvalues <- qPermutvalues <- matrix(NA, ncol = NumComps, nrow = nrow(Data), dimnames = list(rows = rownames(Data), cols = paste("Permutation p-values", 1:(NumComps))))
    lratios <- NULL
  cat("Running rank products and permutations tests ...\n")
  pb <- txtProgressBar(0.9, NumComps)

  for (vs in 1:NumComps) {
    if (!is.null(getDefaultReactiveDomain())) {
      setProgress(0.1 + 0.3 / (NumComps) * vs, detail = paste("tests for comparison", vs, "of", NumComps))
    }
    tData <- Data[, Reps == RRCateg[1, vs]]
    trefData <- Data[, Reps == RRCateg[2, vs]]

    ## t-test
    ttest_out <- ttest_unpaired(tData, trefData)
    ptvalues[, vs ] <- ttest_out$ptvalues
    qtvalues[, vs ] <- ttest_out$qtvalues


    ## rank products
    rp_out <- rp_unpaired(tData, trefData, NumReps)
    pRPvalues[, vs ] <- rp_out$pRPvalues
    qRPvalues[, vs ] <- rp_out$qRPvalues

    ## Permutation tests
    perm_out <- perm_unpaired(tData, trefData, NumReps)
    pPermutvalues[, vs ] <- perm_out$pPermutvalues
    qPermutvalues[, vs ] <- perm_out$qPermutvalues

    lratios <- cbind(lratios, rowMeans(Data[, Reps == RRCateg[1, vs]], na.rm = T) - rowMeans(Data[, Reps == RRCateg[2, vs]], na.rm = T))
    setTxtProgressBar(pb, vs)
  }
  close(pb)

  return(list(
    lratios = lratios, ptvalues = ptvalues, plvalues = plvalues, pRPvalues = pRPvalues, pPermutvalues = pPermutvalues,
    qtvalues = qtvalues, qlvalues = qlvalues, qRPvalues = qRPvalues, qPermutvalues = qPermutvalues, Sds = Sds))
}


# Separate function of calculation of Miss tests. This happens between pairs of tests but does not
# have distinction for pairwise testing
#
# Calculates the missingness statistics for a given dataset and reference matrix.
#
# Args:
#   Data: The dataset to calculate missingness statistics for.
#   RR: The reference matrix specifying the pairs of tests to compare.
#   NumCond: The number of conditions in the dataset.
#   NumReps: The number of replicates per condition.
#
# Returns:
#   A list containing the calculated p-values and q-values for missingness statistics.
MissingStatsDesign <- function(Data, RR, NumCond, NumReps) {
  Reps <- rep(1:NumCond, NumReps)

  NumComps <- ncol(RR) / NumReps
  RRCateg <- RR[, 1:NumComps, drop = F]


  pNAvalues <- matrix(NA, ncol = NumComps, nrow = nrow(Data), dimnames = list(rows = rownames(Data), cols = 1:(NumComps)))
  qNAvalues <- matrix(NA, ncol = NumComps, nrow = nrow(Data), dimnames = list(rows = rownames(Data), cols = 1:(NumComps)))
  cat("Running Miss test ...\n")
  pb <- txtProgressBar(0.9, NumComps)

  for (vs in 1:NumComps) {
    tData <- Data[, Reps == RRCateg[2, vs]]
    trefData <- Data[, Reps == RRCateg[1, vs]]
    tCompDat <- cbind(tData, trefData)
    qs <- quantile(tCompDat, probs = seq(0, 1, 0.01), na.rm = T)
    pvals <- statis <- matrix(NA, nrow(tCompDat), ncol = length(qs))
    for (q in qs) {
      tCompDat[tCompDat < q] <- NA
      NAPDistr <- MissValPDistr(NumReps, sum(is.na(tCompDat)) / (nrow(tCompDat) * 2 * NumReps))
      statis[, which(q == qs)] <- (rowSums(!is.na(tCompDat[, 1:NumReps])) - rowSums(!is.na(tCompDat[, (NumReps + 1):(2 * NumReps)])))
      pvals[, which(q == qs)] <- NAPDistr[abs(statis[, which(q == qs)]) + 1]
    }
    pNAvalues[, vs] <- rowMins(pvals) * (NumReps + 1)
    # qNAvalues[,vs-1] <- qvalue(pNAvalues[,vs-1],lambda=seq(0.1,max(pNAvalues[,vs-1]),length=100))$qvalue
    qNAvalues[, vs] <- p.adjust(pNAvalues[, vs], method = "BH")
    setTxtProgressBar(pb, vs)
  }
  close(pb)
  # print(head(pNAvalues))
  pNAvalues[pNAvalues > 1] <- 1

  return(list(pNAvalues = pNAvalues, qNAvalues = qNAvalues))
}


## Function: MissingStats
## Description: Calculates the p-values and q-values for missing values in a dataset.
## Parameters:
##   - Data: The input dataset.
##   - NumCond: The number of conditions in the dataset.
##   - NumReps: The number of replicates per condition.
## Returns: A list containing the p-values and q-values for missing values.
MissingStats <- function(Data, NumCond, NumReps) {
  Reps <- rep(1:NumCond, NumReps)
  pNAvalues <- matrix(NA, ncol = NumCond - 1, nrow = nrow(Data), dimnames = list(rows = rownames(Data), cols = 1:(NumCond - 1)))
  qNAvalues <- matrix(NA, ncol = NumCond - 1, nrow = nrow(Data), dimnames = list(rows = rownames(Data), cols = 1:(NumCond - 1)))
  for (vs in 2:NumCond) {
    tData <- Data[, Reps == vs]
    trefData <- Data[, Reps == 1]
    tCompDat <- cbind(tData, trefData)
    qs <- quantile(tCompDat, probs = seq(0, 1, 0.01), na.rm = T)
    pvals <- statis <- matrix(NA, nrow(tCompDat), ncol = length(qs))
    for (q in qs) {
      tCompDat[tCompDat < q] <- NA
      NAPDistr <- MissValPDistr(NumReps, sum(is.na(tCompDat)) / (nrow(tCompDat) * 2 * NumReps))
      statis[, which(q == qs)] <- (rowSums(!is.na(tCompDat[, 1:NumReps])) - rowSums(!is.na(tCompDat[, (NumReps + 1):(2 * NumReps)])))
      pvals[, which(q == qs)] <- NAPDistr[abs(statis[, which(q == qs)]) + 1]
    }
    pNAvalues[, vs - 1] <- rowMins(pvals) * (NumReps + 1)
    qNAvalues[, vs - 1] <- p.adjust(pNAvalues[, vs - 1], method = "BH")
  }

  pNAvalues[pNAvalues > 1] <- 1

  return(list(pNAvalues = pNAvalues, qNAvalues = qNAvalues))
}

# Function to determine "optimal" fold-change and q-value thresholds
# fdrtool and use hc.threshold and fc-threshold from ratios (1 standard deviation)
#
# Parameters:
#   - Pvalue: A matrix of p-values for each condition
#   - LogRatios: A matrix of log ratios for each condition
#
# Returns:
#   - A vector containing the mean fold-change threshold and mean q-value threshold
FindFCandQlimAlternative <- function(Pvalue, LogRatios) {
  BestComb <- c(0, 0)
  NumCond <- ncol(LogRatios) + 1
  Pvalue[is.na(Pvalue)] <- 1

  BestHCs <- BestFCs <- vector(, ncol(LogRatios))

  for (i in 1:ncol(LogRatios)) {
    pvals <- Pvalue[, (i - 1) * 5 + 1]
    BestHCs[i] <- hc.thresh(pvals[pvals < 1], plot = F)
    BestFCs[i] <- sd(LogRatios[, i], na.rm = T)
  }

  print(BestHCs)
  print(BestFCs)

  # Calculate mean of all estimated thresholds

  return(c(mean(BestFCs), mean(BestHCs[i])))
}

# Function to determine "optimal" fold-change and q-value thresholds
#
# This function takes in two parameters: Qvalue and LogRatios. Qvalue is a matrix of q-values, and LogRatios is a matrix of log ratios.
# The function calculates the "optimal" fold-change and q-value thresholds by maximizing the percental output of features commonly found for limma, rank products, permutation, and NA tests.
#
# The function first sets any NA values in Qvalue to 1. It then determines the smallest q-value and selects a range of q-values around it.
# Next, the function runs over different fold-change (FC) thresholds and performs tests for each threshold. For each FC threshold, it modifies the Qvalue matrix by setting values to 1 if the corresponding log ratio falls within the FC threshold.
# The function then runs over the range of q-values and calculates the distribution of significant features for each q-value threshold. It calculates the mean of the distributions and keeps track of the combination of FC and q-value thresholds that yield the highest mean.
# Finally, the function returns the combination of FC and q-value thresholds that resulted in the highest mean distribution of significant features.
#
# Example usage:
# Qvalue <- matrix(c(0.01, 0.02, 0.03, 0.04, 0.05), nrow = 5, ncol = 3)
# LogRatios <- matrix(c(1.2, 0.8, 1.5, -0.5, 0.2, 0.9, -1.1, 0.7, 1.8, -0.9, 0.3, 1.1), nrow = 4, ncol = 3)
# thresholds <- FindFCandQlim(Qvalue, LogRatios)
# print(thresholds)
# Output: [1] 0.5 0.1
#
# Note: This function requires the 'matrixStats' and 'parallel' packages to be installed.
FindFCandQlim <- function(Qvalue, LogRatios) {
  BestComb <- c(0, 0)
  BestRegs <- 0
  NumCond <- ncol(LogRatios) + 1

  Qvalue[is.na(Qvalue)] <- 1

  smallestq <- signif(min(Qvalue, na.rm = T))
  qrange <- c(0.1, 0.2, 0.5) * 10^(rep(-10:0, each = 3))
  qrange <- qrange[which.min(abs(smallestq - qrange)):(length(qrange) - 2)]

  # Run over different FC thresholds
  fcRange <- seq(0, max(abs(range(LogRatios, na.rm = T))), length = 100)
  cl <- makeCluster(NumThreads)
  clusterExport(cl = cl, varlist = c("Qvalue", "NumCond", "LogRatios", "qrange"), envir = environment())
  clusterEvalQ(cl = cl, library(matrixStats))
  # print(head(Qvalue))
  BestVals <- parallel::parLapply(cl, fcRange, function(fc) {
    # range of tests to consider:
    for (t in c(1, 2, 4)) {
      tvals <- Qvalue[, (NumCond - 1) * t + 1:(NumCond - 1)]
      tvals[LogRatios < fc & LogRatios > -fc] <- 1
      Qvalue[, (NumCond - 1) * t + 1:(NumCond - 1)] <- tvals
    }
    # Run over range of q-values
    for (qlim in qrange) {
      alldistr <- vector("numeric", NumCond - 1)
      for (t in 1:(NumCond - 1)) {
        distr <- table(rowSums(Qvalue[, seq(t, ncol(Qvalue), NumCond - 1)] < qlim, na.rm = T))
        allregs <- sum(distr[2:length(distr)], na.rm = T)
        # print(distr)
        if (length(distr) > 1 & !is.na(distr["4"])) {
          alldistr[t] <- distr["4"] / allregs
        }
      }
      # print(alldistr)
      if (mean(alldistr) > BestRegs) {
        BestRegs <- mean(alldistr)
        BestComb <- c(fc, qlim)
        # print(BestRegs)
      }
    }
    c(BestRegs, BestComb)
  })
  stopCluster(cl)

  BestRegs <- 0
  for (i in 1:length(BestVals)) {
    if (BestRegs < BestVals[[i]][1]) {
      BestComb <- BestVals[[i]][2:3]
      BestRegs <- BestVals[[i]][1]
    }
    # print(BestVals)
  }
  return(BestComb)
}


# calculated common q-value over different tests. Hommel methods gives
# upper bound for p-values coming from independent or positively dependent tests
#
# Parameters:
#   Qvalue: A matrix of q-values for multiple tests
#   NumComps: The number of comparisons being made
#   NumTests: The number of tests being performed
#
# Returns:
#   A matrix of unified q-values calculated using the Hommel method
UnifyQvals <- function(Qvalue, NumComps, NumTests) {
  cat("Calculating PolySTest FDRs ...\n")
  UnifiedQvalue <- matrix(NA, ncol = NumComps, nrow = nrow(Qvalue))
  for (i in 1:(NumComps)) {
    UnifiedQvalue[, i] <- colMins(apply(Qvalue[, seq(i, ncol(Qvalue) - NumComps, NumComps)], 1, p.adjust, "hommel"), na.rm = T)
  }
  UnifiedQvalue
}

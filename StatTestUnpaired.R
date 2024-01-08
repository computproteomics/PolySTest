#' Unpaired
#'
#' This function performs significance analysis for unpaired data.
#'
#' @param Data A matrix of gene expression data.
#' @param NumCond The number of conditions.
#' @param NumReps The number of replicates per condition.
#'
#' @return A list containing various results of the significance analysis, including lratios, 
#' ptvalues, plvalues, pRPvalues, pPermutvalues, qtvalues, qlvalues, qRPvalues, qPermutvalues, and Sds.
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
#' @details The function adds columns from the randomized full set to reach the minimum number of
#' permutation columns (NumPermCols) replicates. It randomizes the sign as well to avoid 
#' tendencies to one or the other side. In the unpaired case, it also normalizes by the 
#' mean of the entire sample to avoid strange effects. The function then performs permutation 
#' testing using parallel computing, and calculates the p-values and q-values based on the permutation results.
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

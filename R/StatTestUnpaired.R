#' PolySTest for unpaired tests
#'
#' Combining the power of different statistical tests
#'
#' @param fulldata A SummarizedExperiment or derived object that contains
#' the quantitative data as required for PolySTest
#' @param allComps A matrix containing the reference matrix specifying the
#' pairs of conditions to compare (each comparison given as separate row)
#'
#' @return SummarizedExperiment with added columns for p-values and q-values
#' in rowData
#'
#' @examples
#' # Assuming 'fulldata' is a SummarizedExperiment object
#' # and 'allComps' is a matrix containing the reference matrix
#' # specifying the pairs of conditions to compare
#' PolySTest_unpaired(fulldata, allComps)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
PolySTest_unpaired <- function(fulldata, allComps) {

  # check fulldata
  check_for_polystest(fulldata)

  # Extracting the assay data
  Data <- assay(fulldata, "quant")

  # Extracting the metadata
  NumReps <- metadata(fulldata)$NumReps
  NumCond <- metadata(fulldata)$NumCond

  # Extrating contrast details
  Reps <- rep(1:NumCond, NumReps)
  conditions <- unique(colData(fulldata)$Condition)
  NumComps <- nrow(allComps)
  # Make indices for pairings
  RRCateg <- matrix(NA, ncol = nrow(allComps), nrow = 2)
  for (i in 1:nrow(allComps)) {
    RRCateg[1, i] <- as.numeric(which(conditions == allComps[i, 1]))
    RRCateg[2, i] <- as.numeric(which(conditions == allComps[i, 2]))
  }

  # Normalize row-wise by mean
  Data <- Data - rowMeans(Data, na.rm = T)

  # Prepare output data
  tests <-  c("limma", "Miss_Test", "t-test", "rank_products", "permutation_test")
  p_values <- q_values <- matrix(NA, nrow=nrow(Data), ncol=length(tests)*NumComps)
  rownames(p_values) <- rownames(q_values) <- rownames(Data)
  colnames(p_values) <- paste0("p-values_", rep(tests, each=NumComps), "_", rep(1:NumComps, length(tests)))
  colnames(q_values) <- paste0("q-values_", rep(tests, each=NumComps), "_", rep(1:NumComps, length(tests)))

  ## limma
  cat("Running limma tests\n")
  lm_out <- limma_unpaired(Data, NumCond, NumReps, RRCateg)
  p_values[, grep("p-values_limma", colnames(p_values))] <- lm_out$plvalues
  q_values[, grep("q-values_limma", colnames(q_values))] <- lm_out$qlvalues
  Sds <- lm_out$Sds
  cat("limma completed\n")

  cat("Running Miss test\n")
  MissingStats <- MissingStatsDesign(Data, RRCateg, NumCond, NumReps)
  p_values[, grep("p-values_Miss_Test", colnames(p_values))] <- MissingStats$pNAvalues
  q_values[, grep("q-values_Miss_Test", colnames(q_values))] <- MissingStats$qNAvalues
  cat("Miss test completed\n")

  ## rank products + t-test
  lratios <- NULL
  cat("Running rank products and permutations tests ...\n")
  pb <- txtProgressBar(0.9, NumComps)

  for (vs in 1:NumComps) {
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::setProgress(0.1 + 0.3 / (NumComps) * vs, detail = paste("tests for comparison", vs, "of", NumComps))
    }
    tData <- Data[, Reps == RRCateg[1,  vs]]
    trefData <- Data[, Reps == RRCateg[2, vs]]

    ## t-test
    ttest_out <- ttest_unpaired(tData, trefData)
    p_values[, grep("p-values_t-test", colnames(p_values))[vs]] <- ttest_out$ptvalues
    q_values[, grep("q-values_t-test", colnames(q_values))[vs]] <- ttest_out$qtvalues

    ## rank products
    rp_out <- rp_unpaired(tData, trefData)
    p_values[, grep("p-values_rank_products",  colnames(p_values))[vs]] <- rp_out$pRPvalues
    q_values[, grep("q-values_rank_products", colnames(q_values))[vs]] <- rp_out$qRPvalues

    ## Permutation tests
    perm_out <- perm_unpaired(tData, trefData)
    p_values[, grep("p-values_permutation_test", colnames(p_values))[vs]] <- perm_out$pPermutvalues
    q_values[, grep("q-values_permutation_test", colnames(q_values))[vs]] <- perm_out$qPermutvalues

    lratios <- cbind(lratios, rowMeans(Data[, Reps == RRCateg[1, vs]], na.rm = T) - rowMeans(Data[, Reps == RRCateg[2, vs]], na.rm = T))
    setTxtProgressBar(pb, vs)
  }
  close(pb)

  # Prepare output data
  fulldata <- prepare_output_data(fulldata, p_values, q_values, lratios, tests, allComps)

return(fulldata)
}

#' Perform unpaired limma analysis
#' This function performs unpaired limma analysis on Data.
#' @param Data A matrix of gene expression data.
#' @param NumCond The number of conditions in the experiment.
#' @param NumReps The number of replicates per condition.
#' @param RRCateg A matrix specifying the conditons to be compared
#' @return A list containing the following results:
#'  - plvalues: The p-values from limma tests.
#'  - qlvalues: The q-values from limma tests.
#'  - Sds: The standard deviations of the Bayesian linear model.
#'  @details This function performs unpaired limma analysis on Data. It calculates the p-values and q-values for each row, indicating the significance of the difference between the two datasets.
#' @keywords limma unpaired analysis
#' @export
#' @import limma
#' @import qvalue
#' @examples
#'  dataMatrix <- matrix(rnorm(100), nrow = 10)
#' colData <- DataFrame(Condition = rep(c("A", "B"), each = 5))
#' rowData <- DataFrame(Gene = paste("Gene", 1:10))
#' fulldata <- SummarizedExperiment(assays = list(quant = dataMatrix),
#'                                  colData = colData, rowData = rowData)
#' metadata(fulldata) <- list(NumCond = 2, NumReps = 5)
#' #  Specifying comparisons
#' allComps <- matrix(c("A", "B"), ncol = 2, byrow = TRUE)
#' # Run function
#' results <- PolySTest_unpaired(fulldata, allComps)
#' # found differentially regulated features
#'
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
  contrast.matrix <- limma::makeContrasts(contrasts = contrasts, levels = design)
  lm.fitted <- limma::lmFit(Data, design)

  lm.contr <- limma::contrasts.fit(lm.fitted, contrast.matrix)
  lm.bayes <- limma::eBayes(lm.contr)
  topTable(lm.bayes)
  plvalues <- lm.bayes$p.value
  qlvalues <- matrix(NA, nrow = nrow(plvalues), ncol = ncol(plvalues), dimnames = dimnames(plvalues))
  # qvalue correction
  for (i in 1:ncol(plvalues)) {
    tqs <- qvalue::qvalue(na.omit(plvalues[, i]))$qvalues
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
#' @import qvalue
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
  tqs <- qvalue::qvalue(na.omit(ptvalues))$qvalues
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
#' rp_unpaired(tData, trefData)
#'
#' @export
#' @import matrixStats
#' @import parallel
#'
rp_unpaired <- function(tData, trefData) {
  # calculate NumRPPairs random pairing combinations and then take mean of p-values
  if (exists("NumRPPairs")) {
    NumRPPairs <- NumRPPairs
  } else {
    NumRPPairs <- 100
  }

  NumReps <- ncol(tData)
  tpRPvalues <- matrix(NA,
                       ncol = NumRPPairs, nrow = nrow(tData),
                       dimnames = list(rows = rownames(tData), cols = 1:NumRPPairs)
  )
  NumThreads <- get_numthreads()
  cl <- parallel::makeCluster(NumThreads)
  parallel::clusterExport(cl = cl, varlist = c("NumReps", "tData", "trefData", "RPStats"), envir = environment())
  parallel::clusterEvalQ(cl = cl, library(matrixStats))

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
#' result <- perm_unpaired(tData, trefData)
#'
#' @export
#' @import matrixStats
#' @import parallel
#' @import qvalue
perm_unpaired <- function(tData, trefData) {
  NumReps <- ncol(tData)
  # if there is an object NumPermCols, use it, otherwise use default value
  if (exists("NumPermCols")) {
    NumPermCols <- NumPermCols
  } else {
    NumPermCols <- 7
  }
  # if there is an object NumTests, use it, otherwise use default value
  if (exists("NumTests")) {
    NTests <- NumTests
  } else {
    NTests <- 1000
  }

  # add columns from randomized full set to reach min. NumPermCols replicates
  # randomizing also sign to avoid tendencies to one or the other side
  # In the unpaired case, also normalize by mean of the entire sample to avoid strange effects
  tData <- tData - mean(as.numeric(unlist(tData)), na.rm = T)
  trefData <- trefData - mean(as.numeric(unlist(trefData)), na.rm = T)
  if (ncol(tData) * 2 < NumPermCols) {
    AddDat <- matrix(sample(as.vector(unlist(tData)),
                            (NumPermCols - ncol(tData)) * nrow(tData), replace = T),
                     nrow = nrow(tData))
    PermData <- cbind(tData, AddDat)
    AddDat <- matrix(sample(as.vector(unlist(trefData)), (NumPermCols - ncol(trefData)) * nrow(trefData), replace = T), nrow = nrow(trefData))
    PermFullData <- cbind(PermData, trefData, AddDat)
  } else {
    PermFullData <- cbind(tData, trefData)
  }
  RealStats <- StatsForPermutTest(as.matrix(cbind(trefData, tData)), Paired = F)
  # print(head(PermFullData))

  NumThreads <- get_numthreads()
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
  tqs <- qvalue::qvalue(na.omit(pPermutvalues))$qvalues
  qPermutvalues[names(tqs)] <- tqs

  return(list(pPermutvalues = pPermutvalues, qPermutvalues = qPermutvalues))
}

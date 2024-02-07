#' PolySTest for paired tests
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
#' PolySTest_paired(fulldata, allComps)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
PolySTest_paired <- function(fulldata, allComps) {
  # check fulldata
  check_for_polystest(fulldata)

  # Extracting the assay data
  dat <- assay(fulldata, "quant")

  # Extracting the metadata
  NumReps <- metadata(fulldata)$NumReps
  NumCond <- metadata(fulldata)$NumCond

  # Create the ratios matrix
  MAData <- create_ratio_matrix(fulldata, allComps)
  MAReps <- rep(1:nrow(allComps), NumReps)
  ncomps <- nrow(allComps)
  # Make indices for pairings
  RRCateg <- matrix(NA, nrow = nrow(allComps), ncol = 2)
  for (i in 1:nrow(allComps)) {
    RRCateg[i, 1] <- as.numeric(which(conditions == allComps[i, 1]))
    RRCateg[i, 2] <- as.numeric(which(conditions == allComps[i, 2]))
  }


  # Prepare output data
  tests <-  c("limma", "Miss_Test", "t-test", "rank_products", "permutation_test")
  p_values <- q_values <- matrix(NA, nrow=nrow(MAData), ncol=length(tests)*ncomps)
  rownames(p_values) <- rownames(q_values) <- rownames(MAData)
  colnames(p_values) <- paste0("p-values_", rep(tests, each=ncomps), "_", rep(1:ncomps, length(tests)))
  colnames(q_values) <- paste0("q-values_", rep(tests, each=ncomps), "_", rep(1:ncomps, length(tests)))

  cat("Running paired tests\n")
  ## limma with ratios
  limma_out <- limma_paired(MAData, ncomps, NumReps)
  p_values[, grep("p-values_limma", colnames(p_values))] <- limma_out$plvalues
  q_values[, grep("q-values_limma", colnames(q_values))] <- limma_out$qlvalues
  Sds <- limma_out$Sds

  cat("limma completed\n")

  cat("Running Miss test\n")
  MissingStats <- MissingStatsDesign(dat, RRCateg, NumCond, NumReps)
  p_values[, grep("p-values_Miss_Test", colnames(p_values))] <- MissingStats$pNAvalues
  q_values[, grep("q-values_Miss_Test", colnames(q_values))] <- MissingStats$qNAvalues
  cat("Miss test completed\n")

  cat("Running rank products and permutations tests ...\n")
  lratios <- NULL
  pb <- txtProgressBar(0.9, NumCond)
  for (vs in 1:ncomps) {
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::setProgress(0.1 + 0.3 / ncomps * vs, detail = paste("tests for comparison", vs, "of", ncomps))
    }

    tMAData <- MAData[, MAReps == vs]

    ## t-tests
    ttest_out <- ttest_paired(tMAData)
    p_values[, grep("p-values_t-test", colnames(p_values))[vs]] <- ttest_out$ptvalues
    q_values[, grep("q-values_t-test", colnames(q_values))[vs]] <- ttest_out$qtvalues

    ## rank products
    # Up
    RPMAUp_pvalues <- RPStats(tMAData, NumReps)
    # Down
    RPMADown_pvalues <- RPStats(-tMAData, NumReps)
    ttt <- rowMins(cbind(RPMAUp_pvalues, RPMADown_pvalues), na.rm = T) * 2
    ttt[ttt > 1] <- 1
    p_values[names(RPMAUp_pvalues), grep("p-values_rank_products", colnames(p_values))[vs]] <- ttt
    tqs <- p.adjust(na.omit(ttt), method = "BH")
    q_values[names(tqs), grep("q-values_rank_products", colnames(q_values))[vs]] <- tqs

    ## Permutation tests
    perm_out <- permtest_paired(tMAData)
    p_values[, grep("p-values_permutation_test", colnames(p_values))[vs]] <- perm_out$pPermutvalues
    q_values[, grep("q-values_permutation_test", colnames(q_values))[vs]] <- perm_out$qPermutvalues

    lratios <- cbind(lratios, rowMeans(MAData[, MAReps == i], na.rm = T))
    setTxtProgressBar(pb, vs)
  }
  cat("rank products and permutation test completed\n")
  close(pb)

  # Prepare output data
  fulldata <- prepare_output_data(fulldata, p_values, q_values, lratios, tests, allComps)

  return(fulldata)
}

#' Perform paired limma analysis
#'
#' This function performs paired limma analysis on MAData.
#'
#' @param MAData A ratio matrix of gene expression data. The rows are genes and the columns are samples.
#' Replicates must be grouped together.
#' @param NumCond The number of ratios to check vs zero level
#' @param NumReps The number of replicates per condition.
#'
#' @return A list containing the p-values and q-values.
#' @keywords limma paired analysis
#' @import limma
#' @import qvalue
#' @export
limma_paired <- function(MAData, NumCond, NumReps) {
  MAReps <- rep(1:NumCond, NumReps)
  ## limma with ratios
  design <- plvalues <- NULL
  for (c in (1:(NumCond))) {
    design <- cbind(design, as.numeric(MAReps == c))
  }
  lm.fittedMA <- limma::lmFit(MAData, design)
  lm.bayesMA <- limma::eBayes(lm.fittedMA)
  topTable(lm.bayesMA)
  plvalues <- lm.bayesMA$p.value
  qlvalues <- matrix(NA, nrow = nrow(plvalues), ncol = ncol(plvalues), dimnames = dimnames(plvalues))
  # qvalue correction
  for (i in 1:ncol(plvalues)) {
    tqs <- qvalue::qvalue(na.omit(plvalues[, i]))$qvalues
    qlvalues[names(tqs), i] <- tqs
  }
  return(list(plvalues = plvalues, qlvalues = qlvalues, Sds = sqrt(lm.bayesMA$s2.post)))
}

#' Perform paired t-tests
#'
#' This function performs row-wise paired t-tests on a data frame containing
#' log-fold changes
#'
#' @param tMAData A matrix of data for running row-wise t-tests
#'
#' @return A list containing the p-values and q-values (qvalue package)
#' @keywords t-test paired analysis
#' @export
#' @import qvalue
#' @examples
#' tMAData <- matrix(rnorm(100), nrow = 10)
#' tout <- ttest_paired(tMAData)
#' head(tout$qtvalues)
ttest_paired <- function(tMAData) {
  ## t-tests
  ptvalues <- sapply(1:nrow(tMAData), function(pep) {
    ifelse(sum(!is.na(tMAData[pep, ])) > 1,
           t.test(tMAData[pep, ])$p.value,
           NA
    )
  })
  names(ptvalues) <- rownames(tMAData)
  # Storey FDR correction
  tqs <- qvalue::qvalue(na.omit(ptvalues))$qvalues
  qtvalues <- rep(NA, length(ptvalues))
  names(qtvalues) <- names(ptvalues)
  qtvalues[names(tqs)] <- tqs

  return(list(ptvalues = ptvalues, qtvalues = qtvalues))
}

#' Perform permutation tests
#'
#' This function performs permutation tests on the given tMAData
#' The permutation tests determine an empirical null distribution of t-values for
#' the p-value calculation
#'
#' @param tMAData A matrix of data for running permutation tests
#'
#' @return A list containing the p-values and q-values (Benjamini-Hochberg)
#' @keywords permutation paired analysis
#' @export
#' @import parallel
#' @examples
#' tMAData <- matrix(rnorm(100), nrow = 10)
#' tout <- permtest_paired(tMAData)
#' head(tout$qPermutvalues)
permtest_paired <- function(tMAData) {
  ## Permutation tests
  NumReps <- ncol(tMAData)
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

  # if necessary, add columns from randomized full set to reach min. NumPermCols replicates
  # randomizing also sign to avoid tendencies to one or the other side
  if (ncol(tMAData) < NumPermCols) {
    AddDat <- matrix(sample(as.vector(tMAData), (NumPermCols - ncol(tMAData)) * nrow(tMAData), replace = T), nrow = nrow(tMAData))
    PermMAData <- cbind(tMAData, AddDat)
  } else {
    PermMAData <- tMAData
  }
  # calculate t-values for real data
  RealStats <- StatsForPermutTest(tMAData, Paired = T)
  # run in parallel to speed up
  NumThreads <- get_numthreads()
  cl <- parallel::makeCluster(NumThreads)
  parallel::clusterExport(cl = cl, varlist = c("NumReps", "PermMAData", "RPStats", "StatsForPermutTest"), envir = environment())
  parallel::clusterEvalQ(cl = cl, library(matrixStats))

  # calculate t-values for permuted data
  PermutOut <- parallel::parLapply(cl, seq_len(NTests), function(x) {
    indat <- apply(
      PermMAData, 1,
      function(y) sample(y, NumReps) * sample(c(1, -1), NumReps, replace = T)
    )
    StatsForPermutTest(t(indat), T)
  })
  parallel::stopCluster(cl)

  # calculate p-values
  PermutOut <- matrix(unlist(PermutOut), nrow = nrow(tMAData))
  pPermutvalues <- apply(cbind(RealStats, PermutOut), 1, function(x) ifelse(is.na(x[1]), NA, (1 + sum(x[1] < x[-1], na.rm = T)) / (sum(!is.na(x)))))
  qPermutvalues <- rep(NA, length(pPermutvalues))
  names(qPermutvalues) <- names(pPermutvalues)
  # Benjamini-Hochberg FDR correction
  tqs <- p.adjust(na.omit(pPermutvalues), method = "BH")
  qPermutvalues[names(tqs)] <- tqs

  return(list(pPermutvalues = pPermutvalues, qPermutvalues = qPermutvalues))
}

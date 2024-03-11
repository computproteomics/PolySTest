#' PolySTest for Paired Tests
#'
#' Integrates various statistical tests for analyzing paired data within
#' a `SummarizedExperiment` framework. It compares pairs of conditions
#' specified by the user to calculate p-values and q-values for each comparison.
#'
#' @param fulldata A `SummarizedExperiment` or derived object containing
#'                 the quantitative data required for PolySTest analysis.
#' @param allComps A matrix specifying pairs of conditions to compare.
#'                 Each row represents a pair for comparison.
#' @param statTests A character vector specifying the statistical tests to
#'                  be applied. Available tests include "limma", "Miss_Test",
#'                  "t-test", "rank_products", and "permutation_test". The function
#'                  will perform each specified test and integrate the results.
#'
#' @details Executes specified statistical tests on the dataset contained in
#'          `fulldata` using the condition pairs outlined in `allComps`. Calculates
#'          p-values and q-values for each genomic feature (e.g., genes, proteins)
#'          included in the analysis. Results are added to the `rowData` of the
#'          `SummarizedExperiment` object, enhancing it with detailed statistics
#'          from the paired tests.
#'
#' @return A `SummarizedExperiment` object augmented with p-values and q-values
#'         in its `rowData`, reflecting the outcomes of the specified statistical analyses.
#'
#' @examples
#' library(SummarizedExperiment)
#'
#' # Mock quantitative data and metadata for samples
#' quantData <- matrix(rnorm(2000), nrow=200, ncol=10)
#' colnames(quantData) <- c(paste("Sample", 1:5, "_Condition_A", sep=""),
#'                          paste("Sample", 1:5, "_Condition_B", sep=""))
#' rownames(quantData) <- paste("Gene", 1:200)
#' sampleMetadata <- data.frame(Condition=rep(c("A", "B"), each=5))
#'
#' # Creating the SummarizedExperiment object
#' fulldata <- SummarizedExperiment(assays=list(quant=quantData),
#'                                  colData=sampleMetadata)
#' metadata(fulldata) <- list(NumReps=5, NumCond=2)
#'
#' # Specifying pairs of conditions to compare
#' allComps <- matrix(c("A", "B"), ncol=2, byrow=TRUE)
#'
#' # Specify statistical tests to apply
#' statTests <- c("limma", "t_test", "rank_products")
#'
#' # Running PolySTest for paired comparisons
#' results <- PolySTest_paired(fulldata, allComps, statTests)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
PolySTest_paired <- function(fulldata, allComps, statTests = c("limma", "Miss_Test", "t_test", "rank_products", "permutation_test")) {

  # Check if the fulldata object meets the requirements
  check_for_polystest(fulldata)

  # Extracting the assay data
  Data <- assay(fulldata, "quant")

  # Extract metadata from the SummarizedExperiment object
  NumReps <- metadata(fulldata)$NumReps
  NumCond <- metadata(fulldata)$NumCond


  # Validate specified statistical tests
  if (!all(statTests %in% c("limma", "Miss_Test", "t_test", "rank_products", "permutation_test"))) {
    stop("Invalid statistical test specified. Available tests are: limma, Miss_Test, t_test, rank_products, permutation_test")
  }
  tests <- statTests


  # Create the ratios matrix
  MAData <- create_ratio_matrix(fulldata, allComps)
  MAReps <- rep(seq_len(nrow(allComps)), NumReps)

  # Extrating contrast details
  Reps <- rep(seq_len(NumCond), NumReps)
  conditions <- unique(colData(fulldata)$Condition)
  NumComps <- nrow(allComps)
  # Make indices for pairings
  RRCateg <- matrix(NA, ncol = nrow(allComps), nrow = 2)
  for (i in seq_len(nrow(allComps))) {
    RRCateg[1, i] <- as.numeric(which(conditions == allComps[i, 1]))
    RRCateg[2, i] <- as.numeric(which(conditions == allComps[i, 2]))
  }

  # Prepare output data
  p_values <- q_values <- matrix(NA, nrow=nrow(MAData), ncol=length(tests)*NumComps)
  rownames(p_values) <- rownames(q_values) <- rownames(MAData)
  colnames(p_values) <- paste0("p_values_", rep(tests, each=NumComps), "_", rep(seq_len(NumComps), length(tests)))
  colnames(q_values) <- paste0("q_values_", rep(tests, each=NumComps), "_", rep(seq_len(NumComps), length(tests)))

  Sds <- NULL
  if (any("limma" %in% tests)) {
    message("Running limma tests")
    ## limma with ratios
    limma_out <- limma_paired(MAData, NumComps, NumReps)
    p_values[, grep("p_values_limma", colnames(p_values))] <- limma_out$plvalues
    q_values[, grep("q_values_limma", colnames(q_values))] <- limma_out$qlvalues
    Sds <- limma_out$Sds
    message("limma completed")
  }

  if (any("Miss_Test" %in% tests)) {
    message("Running Miss test")
    MissingStats <- MissingStatsDesign(Data, RRCateg, NumCond, NumReps)
    p_values[, grep("p_values_Miss_Test", colnames(p_values))] <- MissingStats$pNAvalues
    q_values[, grep("q_values_Miss_Test", colnames(q_values))] <- MissingStats$qNAvalues
    message("Miss test completed")
  }

  message("Running rank products and permutations tests ...")
  lratios <- NULL
  if (any("rank_products" %in% tests)) {
    message("Running rank products ...")
  }
  if (any("t_test" %in% tests)) {
    message("Running t-tests ...")
  }
  if (any("permutation_test" %in% tests)) {
    message("Running permutation tests ...")
  }
  pb <- txtProgressBar(0.9, NumCond)
  for (vs in seq_len(NumComps)) {
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::setProgress(0.1 + 0.3 / NumComps * vs, detail = paste("tests for comparison", vs, "of", NumComps))
    }

    tMAData <- MAData[, MAReps == vs]

    ## t-tests
    if (any("t_test" %in% tests)) {

      ttest_out <- ttest_paired(tMAData)
      p_values[, grep("p_values_t_test", colnames(p_values))[vs]] <- ttest_out$ptvalues
      q_values[, grep("q_values_t_test", colnames(q_values))[vs]] <- ttest_out$qtvalues
    }

    ## rank products
    if (any("rank_products" %in% tests)) {
      # Up
      RPMAUp_pvalues <- RPStats(tMAData, NumReps)
      # Down
      RPMADown_pvalues <- RPStats(-tMAData, NumReps)
      ttt <- rowMins(cbind(RPMAUp_pvalues, RPMADown_pvalues), na.rm = T) * 2
      ttt[ttt > 1] <- 1
      p_values[names(RPMAUp_pvalues), grep("p_values_rank_products", colnames(p_values))[vs]] <- ttt
      tqs <- p.adjust(na.omit(ttt), method = "BH")
      q_values[names(tqs), grep("q_values_rank_products", colnames(q_values))[vs]] <- tqs
    }

    ## Permutation tests
    if (any("permutation_test" %in% tests)) {
      perm_out <- permtest_paired(tMAData)
      p_values[, grep("p_values_permutation_test", colnames(p_values))[vs]] <- perm_out$pPermutvalues
      q_values[, grep("q_values_permutation_test", colnames(q_values))[vs]] <- perm_out$qPermutvalues
    }

    lratios <- cbind(lratios, rowMeans(MAData[, MAReps == i], na.rm = T))
    setTxtProgressBar(pb, vs)
  }
  message("rank products and permutation test completed")
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
#'
#' @examples
#' MAData <- matrix(rnorm(600), nrow = 100)
#' NumCond <- 3
#' NumReps <- 2
#' limma_res <- limma_paired(MAData, NumCond, NumReps)
#' head(limma_res$qlvalues)
#'
#' @keywords limma paired analysis
#' @import limma
#' @import qvalue
#' @export
limma_paired <- function(MAData, NumCond, NumReps) {
  MAReps <- rep(seq_len(NumCond), NumReps)
  ## limma with ratios
  design <- plvalues <- NULL
  for (c in (seq_len(NumCond))) {
    design <- cbind(design, as.numeric(MAReps == c))
  }
  lm.fittedMA <- limma::lmFit(MAData, design)
  lm.bayesMA <- limma::eBayes(lm.fittedMA)
  topTable(lm.bayesMA)
  plvalues <- lm.bayesMA$p.value
  qlvalues <- matrix(NA, nrow = nrow(plvalues), ncol = ncol(plvalues), dimnames = dimnames(plvalues))
  # qvalue correction
  for (i in seq_len(ncol(plvalues))) {
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
#' tMAData <- matrix(rnorm(1000), nrow = 100)
#' tout <- ttest_paired(tMAData)
#' head(tout$qtvalues)
ttest_paired <- function(tMAData) {
  ## t-tests
  ptvalues <- sapply(seq_len(nrow(tMAData)), function(pep) {
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

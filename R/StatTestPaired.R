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
#'                  "t-test", "rank_products", and "permutation_test". The
#'                  function
#'                  will perform each specified test and integrate the results.
#'
#' @details Executes specified statistical tests on the dataset contained in
#'          `fulldata` using the condition pairs outlined in `allComps`.
#'          Calculates p-values and q-values for each genomic feature (e.g.,
#'          genes, proteins) included in the analysis. Results are added to the
#'          `rowData` of the
#'          `SummarizedExperiment` object, enhancing it with detailed statistics
#'          from the paired tests.
#'
#' @return A `SummarizedExperiment` object augmented with p-values and q-values
#'         in its `rowData`, reflecting the outcomes of the specified
#'         statistical analyses.
#'
#' @examples
#' library(SummarizedExperiment)
#'
#' # Mock quantitative data and metadata for samples
#' quantData <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#' colnames(quantData) <- c(
#'     paste("Sample", 1:5, "_Condition_A", sep = ""),
#'     paste("Sample", 1:5, "_Condition_B", sep = "")
#' )
#' rownames(quantData) <- paste("Gene", 1:200)
#' sampleMetadata <- data.frame(Condition = rep(c("A", "B"), each = 5))
#'
#' # Creating the SummarizedExperiment object
#' fulldata <- SummarizedExperiment(
#'     assays = list(quant = quantData),
#'     colData = sampleMetadata
#' )
#' metadata(fulldata) <- list(NumReps = 5, NumCond = 2)
#'
#' # Specifying pairs of conditions to compare
#' allComps <- matrix(c("A", "B"), ncol = 2, byrow = TRUE)
#'
#' # Specify statistical tests to apply
#' statTests <- c("limma", "t_test", "rank_products")
#'
#' # Running PolySTest for paired comparisons
#' results <- PolySTest_paired(fulldata, allComps, statTests)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#' @importFrom matrixStats rowMins
PolySTest_paired <- function(fulldata,
                             allComps,
                             statTests = c(
                                 "limma", "Miss_Test", "t_test",
                                 "rank_products",
                                 "permutation_test"
                             )) {
    # Check if the fulldata object meets the requirements
    check_for_polystest(fulldata)
    
    # Extracting the assay data
    Data <- assay(fulldata, "quant")
    
    # Extract metadata from the SummarizedExperiment object
    NumReps <- metadata(fulldata)$NumReps
    NumCond <- metadata(fulldata)$NumCond
    
    
    # Validate specified statistical tests
    if (!all(statTests %in% c(
        "limma", "Miss_Test", "t_test", "rank_products",
        "permutation_test"
    ))) {
        stop("Invalid statistical test specified. Available tests are:
         limma, Miss_Test, t_test, rank_products, permutation_test")
    }
    tests <- statTests
    
    
    # Create the ratios matrix
    MAData <- create_ratio_matrix(fulldata, allComps)
    MAReps <- rep(seq_len(nrow(allComps)), NumReps)
    
    # Extracting contrast details
    conditions <- unique(colData(fulldata)$Condition)
    NumComps <- nrow(allComps)
    # Make indices for pairings
    RRCateg <- matrix(NA, ncol = nrow(allComps), nrow = 2)
    RRCateg[1, ] <- match(allComps[, 1], conditions)
    RRCateg[2, ] <- match(allComps[, 2], conditions)
    
    # Prepare output data
    
    p_values <- q_values <- matrix(NA, nrow = nrow(MAData), 
                                   ncol = length(statTests) * NumComps)
    rownames(p_values) <- rownames(q_values) <- rownames(MAData)
    colnames(p_values) <- paste0("p_values_", rep(statTests, each = NumComps), 
                                 "_", rep(seq_len(NumComps),
                                          length(statTests)))
    colnames(q_values) <- paste0("q_values_", rep(statTests, each = NumComps),
                                 "_", rep(seq_len(NumComps), 
                                          length(statTests)))
    
    
    # Define a list of test functions
    test_funcs <- list(
        limma = function(tMAData) {
            limma_out <- limma_paired(tMAData, NumComps, NumReps)
            list(pvals = limma_out$plvalues, qvals = limma_out$qlvalues, Sds = limma_out$Sds)
        },
        Miss_Test = function(tMAData) {
            MissingStats <- MissingStatsDesign(Data, RRCateg, NumCond, NumReps)
            list(pvals = MissingStats$pNAvalues, qvals = MissingStats$qNAvalues)
        },
        t_test = function(tMAData) {
            ttest_out <- ttest_paired(tMAData)
            list(pvals = ttest_out$ptvalues, qvals = ttest_out$qtvalues)
        },
        rank_products = function(tMAData) {
            RPMAUp_pvalues <- RPStats(tMAData, NumReps)
            RPMADown_pvalues <- RPStats(-tMAData, NumReps)
            ttt <- rowMins(cbind(RPMAUp_pvalues, RPMADown_pvalues), 
                           na.rm = TRUE) * 2
            ttt[ttt > 1] <- 1
            tqs <- rep(NA, length(ttt))
            tqs[!is.na(ttt)] <- p.adjust(na.omit(ttt), method = "BH")
            list(pvals = ttt, qvals = tqs)
        },
        permutation_test = function(tMAData) {
            perm_out <- permtest_paired(tMAData)
            list(pvals = perm_out$pPermutvalues, 
                 qvals = perm_out$qPermutvalues)
        }
    )        
    
    # Run tests that do not depend on comparisons
    Sds <- NULL
    
    for (test in statTests[statTests %in% names(test_funcs)]) {
        if (test %in% c("limma", "Miss_Test")) {
            message("Running", test, "test")
            res <- test_funcs[[test]](MAData)
            p_values[, grep(paste0("p_values_", test), 
                            colnames(p_values))] <- res$pvals
            q_values[, grep(paste0("q_values_", test), 
                            colnames(q_values))] <- res$qvals
            if (test == "limma") Sds <- res$Sds
            message(test, "completed")
        }
    }

    
    lratios <- NULL
    pb <- txtProgressBar(0.9, NumCond)
    for (vs in seq_len(NumComps)) {
        if (!is.null(shiny::getDefaultReactiveDomain())) {
            shiny::setProgress(0.1 + 0.3 / NumComps * vs,
                               detail = paste(
                                   "tests for comparison",
                                   vs, "of", NumComps
                               )
            )
        }
        
        tMAData <- MAData[, MAReps == vs, drop=FALSE]
        
        
        for (test in statTests[statTests %in% c("t_test", "rank_products", "permutation_test")]) {
            res <- test_funcs[[test]](tMAData)
            p_values[, grep(paste0("p_values_", test), 
                            colnames(p_values))[vs]] <- res$pvals
            q_values[, grep(paste0("q_values_", test), 
                            colnames(q_values))[vs]] <- res$qvals
        }        

        
        lratios <- cbind(lratios, rowMeans(MAData[, MAReps == vs], 
                                           na.rm = TRUE))
        
        setTxtProgressBar(pb, vs)
    }
    message("tests completed")
    close(pb)
    
    # Prepare output data
    fulldata <- prepare_output_data(
        fulldata, p_values, q_values,
        lratios, tests, allComps
    )
    
    return(fulldata)
}

#' Perform paired limma analysis
#'
#' This function performs paired limma analysis on MAData.
#'
#' @param MAData A ratio matrix of gene expression data. The rows are genes and
#' the columns are samples.
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
#' @importFrom limma lmFit eBayes topTable
#' @importFrom qvalue qvalue
#' @export
limma_paired <- function(MAData, NumCond, NumReps) {
    MAReps <- rep(seq_len(NumCond), NumReps)
    ## limma with ratios
    design <- plvalues <- NULL
    for (c in (seq_len(NumCond))) {
        design <- cbind(design, as.numeric(MAReps == c))
    }
    lm.fittedMA <- limma::lmFit(MAData, design)
    
    res <- fit_and_getvals(lm.fittedMA)
    
    return(res)    
    
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
#' @importFrom qvalue qvalue
#' @examples
#' tMAData <- matrix(rnorm(1000), nrow = 100)
#' tout <- ttest_paired(tMAData)
#' head(tout$qtvalues)
ttest_paired <- function(tMAData) {
    ## t-tests
    ptvalues <- vapply(seq_len(nrow(tMAData)), function(pep) {
        ifelse(sum(!is.na(tMAData[pep, ])) > 1,
               t.test(tMAData[pep, ])$p.value,
               NA_real_
        )
    }, numeric(1))
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
#' The permutation tests determine an empirical null distribution of t-values
#' for the p-value calculation
#'
#' @param tMAData A matrix of data for running permutation tests
#'
#' @return A list containing the p-values and q-values (Benjamini-Hochberg)
#' @keywords permutation paired analysis
#' @export
#' @importFrom parallel detectCores  makeCluster clusterExport stopCluster parLapply clusterEvalQ
#' 
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
    
    # if necessary, add columns from randomized full set to reach min.
    # NumPermCols replicates randomizing also sign to avoid tendencies to one or
    # the other side
    if (ncol(tMAData) < NumPermCols) {
        AddDat <- matrix(sample(as.vector(tMAData),
                                (NumPermCols - ncol(tMAData)) * nrow(tMAData),
                                replace = TRUE
        ), nrow = nrow(tMAData))
        PermMAData <- cbind(tMAData, AddDat)
    } else {
        PermMAData <- tMAData
    }
    # calculate t-values for real data
    RealStats <- StatsForPermutTest(tMAData, Paired = TRUE)
    # run in parallel to speed up
    NumThreads <- get_numthreads()
    cl <- parallel::makeCluster(NumThreads)
    parallel::clusterExport(
        cl = cl, varlist = c(
            "NumReps", "PermMAData",
            "RPStats", "StatsForPermutTest"
        ),
        envir = environment()
    )
    
    # calculate t-values for permuted data
    PermutOut <- parallel::parLapply(cl, seq_len(NTests), function(x) {
        indat <- apply(
            PermMAData, 1,
            function(y) {
                sample(y, NumReps) * sample(c(1, -1), NumReps,
                                            replace = TRUE
                )
            }
        )
        StatsForPermutTest(t(indat), TRUE)
    })
    parallel::stopCluster(cl)
    
    # calculate p-values
    PermutOut <- matrix(unlist(PermutOut), nrow = nrow(tMAData))
    pPermutvalues <- apply(
        cbind(RealStats, PermutOut), 1,
        function(x) {
            ifelse(is.na(x[1]), NA,
                   (1 + sum(x[1] < x[-1],
                            na.rm = TRUE
                   )) /
                       (sum(!is.na(x)))
            )
        }
    )
    qPermutvalues <- rep(NA, length(pPermutvalues))
    names(qPermutvalues) <- names(pPermutvalues)
    # Benjamini-Hochberg FDR correction
    tqs <- p.adjust(na.omit(pPermutvalues), method = "BH")
    qPermutvalues[names(tqs)] <- tqs
    
    return(list(pPermutvalues = pPermutvalues, qPermutvalues = qPermutvalues))
}

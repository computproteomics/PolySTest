# General settings: Permutations: min 7 replicates, if not available, then
# generate from entire data

#' @importFrom S4Vectors metadata
#' @importFrom grDevices adjustcolor grey.colors rainbow
#' @importFrom graphics abline axis hist layout legend lines
#' @importFrom graphics mtext par text title
#' @importFrom stats model.matrix na.omit p.adjust pgamma quantile runif
#' @importFrom stats sd t.test
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom S4Vectors metadata
NULL

#' Set number of threads (default is 4)
#'
#' This function sets the number of threads for parallel processing. If the
#' environment variable SHINY_THREADS is set, the number of threads is set
#' to the value of SHINY_THREADS.
#'
#' @param threads An integer indicating the number of threads to use.
#' 
#' @return An integer indicating the number of threads to use.
#' 
#' @examples
#' get_numthreads(threads = 4)
#' get_numthreads()
#' @export
#' 
get_numthreads <- function(threads = NULL) {
    NumThreads <- 1
    shiny_threads <- as.numeric(Sys.getenv("SHINY_THREADS"))
    if (!is.na(shiny_threads)) {
        NumThreads <- shiny_threads
        message("Set number of threads to ", NumThreads)
    }
    if (!is.null(threads)) {
        NumThreads <- threads
    }
    return(NumThreads)
}

# Function to select colors for columns
colSelected <- function(col, num, sel, col2) {
    ttt <- rep(col, num)
    ttt[sel] <- col2
    return(ttt)
}


#' Calculate statistics for permutation test
#'
#' This function calculates the statistics for a permutation test based on the
#' input data.
#'
#' @param Data A matrix or data frame containing the data for the test.
#' @param Paired A logical value indicating whether the test is paired or not.
#' @return A numeric vector containing the calculated statistics.
#' @examples
#' Data <- matrix(rnorm(100), ncol = 10)
#' StatsForPermutTest(Data, Paired = FALSE)
#' @importFrom matrixStats rowSds rowVars
#' @export
StatsForPermutTest <- function(Data, Paired) {
    if (Paired) {
        Stats <- rowMeans(Data, na.rm = TRUE) /
            matrixStats::rowSds(Data, na.rm = TRUE)
        Stats <- abs(Stats) * sqrt(rowSums(!is.na(Data)))
    } else {
        NumReps <- ncol(Data) * 0.5
        NumRealReps <- rowSums(!is.na(Data)) * 0.5
        Stats <- rowMeans(Data[, seq_len(NumReps), drop=FALSE], na.rm = TRUE) -
            rowMeans(Data[, seq(NumReps + 1, length.out = NumReps), drop=FALSE],
                     na.rm = TRUE
            ) /
            (sqrt(matrixStats::rowVars(Data[, seq_len(NumReps), drop=FALSE],
                                       na.rm = TRUE
            ) +
                matrixStats::rowVars(
                    Data[, seq(NumReps + 1,
                               length.out = NumReps
                    ), drop=FALSE],
                    na.rm = TRUE
                )))
        Stats <- abs(Stats) * NumRealReps
    }
    return(Stats)
}


#' Calculate the distribution of missing values
#'
#' This function calculates the distribution of missing values for a given
#' number of repetitions and percentage of missing values.
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
    
    # Calculate binomial terms
    it <- seq_len(d+1) - 1
    binTerms <- choose(d, it) * p^(it) * (1 - p)^(rev(it))
    
    # Compute D using vectorized operations
    D <- vapply(it, function(i) sum(binTerms[seq_len(d-i+1)] * binTerms[seq.int(i+1,d+1)]), numeric(1))
    
    # Double the non-zero elements
    D[-1] <- 2 * D[-1]
    
    D
}

#' RPStats Function
#'
#' This function calculates the p-values for the RP (Rank Product) statistic
#' based on the input tRPMAData and the number of replicates (NumReps).
#' The statistcs is one-sided, i.e. it only detects up-regulation.
#'
#' @param tRPMAData A matrix containing the expression data with rows as
#' features and columns as replicates.
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
    
    # This code calculates p-values for RPMADown analysis based on the number of
    # elements (d).
    # It returns a list of p-values for each iteration of d.
    RPMADown_pvalues <- lapply(iterNumEl, function(d) {
        tRPMADown_pvalues <- NULL
        # Subset the RPMAData based on the number of elements (d)
        RPMAData <- tRPMAData[NumElements == d, , drop = FALSE]
        # If d is greater than 1 and the number of columns in RPMAData is
        # greater than the number of columns in tRPMAData
        if (d > 1 && length(as.matrix(RPMAData)) > ncol(tRPMAData)) {
            
            RP.own <- rep(0, nrow(RPMAData))
            names(RP.own) <- rownames(RPMAData)
            
            # Calculate ranks and normalize
            ranks <- apply(RPMAData, 2, rank, na.last = "keep")
            norm_ranks <- apply(ranks, 2, function(x) x / (sum(!is.na(x)) + 1))
            norm_ranks[is.na(norm_ranks)] <- 1
            
            
            # Sum the log of the ranks across columns
            RP.own <- rowSums(log(norm_ranks))
            
            # Calculate RP.own values and p-values
            RP.own <- exp(RP.own)
            tNumFeat <- length(RP.own)
            tRPout <- pgamma(-log(RP.own), d)
            names(tRPout) <- names(RP.own)
            tRPMADown_pvalues <- c(tRPMADown_pvalues, tRPout)
            # Return the calculated p-values
            tRPMADown_pvalues
        } else {
            # If d is less than or equal to 1 or the number of columns in
            # RPMAData is not greater than the number of columns in tRPMAData,
            # return NA values for all rows in RPMAData
            na_out <- as.numeric(rep(NA, nrow(RPMAData)))
            names(na_out) <- rownames(RPMAData)
            return(na_out)
        }
    })
    # reorder to original order
    all_ps <- unlist(RPMADown_pvalues)
    res <- rep(NA, nrow(tRPMAData))
    names(res) <- rownames(tRPMAData)
    res[names(all_ps)] <- all_ps 
    return(res)
}



#' Function of calculation of Miss tests. This happens between full
#' groups and thus does not have distinction for pairwise testing.
#'
#' @param Data A matrix containing the expression data with rows as features and
#' columns as samples
#' @param RRCateg A matrix containing the indices of the conditions to compare
#' (each row is a pairing)
#' @param NumCond The number of conditions in the dataset
#' @param NumReps The number of replicates for each condition (needs to be same
#' for all conditions)
#'
#' @return A list containing the calculated p-values and q-values
#' (Benjamini-Hochbeg) for missingness statistics
#'
#' @examples
#' Data <- matrix(rnorm(120), nrow = 10)
#' # Introduce some missingness
#' Data[sample(seq_len(120), 40)] <- NA
#' RRCateg <- matrix(c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4), nrow = 2, ncol = 6)
#' NumCond <- 4
#' NumReps <- 3
#' res_misstest <- MissingStatsDesign(Data, RRCateg, NumCond, NumReps)
#' head(res_misstest$qNAvalues)
#'
#' @export
MissingStatsDesign <- function(Data, RRCateg, NumCond, NumReps) {
    Reps <- rep(seq_len(NumCond), NumReps)
    
    NumComps <- ncol(RRCateg)
    pNAvalues <- matrix(NA,
                        ncol = NumComps,
                        nrow = nrow(Data),
                        dimnames = list(
                            rows = rownames(Data),
                            cols = seq_len(NumComps)
                        )
    )
    qNAvalues <- matrix(NA,
                        ncol = NumComps, nrow = nrow(Data),
                        dimnames = list(
                            rows = rownames(Data),
                            cols = seq_len(NumComps)
                        )
    )
    
    message("Running Miss test ...")
    pb <- txtProgressBar(0.9, NumComps)
    
    for (vs in seq_len(NumComps)) {
        tData <- Data[, Reps == RRCateg[2, vs], drop=FALSE]
        trefData <- Data[, Reps == RRCateg[1, vs], drop=FALSE]
        tCompDat <- cbind(tData, trefData)
        qs <- quantile(tCompDat, probs = seq(0, 1, 0.01), na.rm = TRUE)
        
        pvals <- vapply(qs, function(q) {
            tCompDatQ <- tCompDat
            tCompDatQ[tCompDatQ < q] <- NA
            NAPDistr <- MissValPDistr(NumReps, sum(is.na(tCompDatQ)) /
                                          (nrow(tCompDatQ) * 2 * NumReps))
            statis <- rowSums(!is.na(tCompDatQ[, seq_len(NumReps)])) -
                rowSums(!is.na(tCompDatQ[, (NumReps + 1):(2 * NumReps)]))
            pvals_q <- NAPDistr[abs(statis) + 1]
            return(pvals_q)
        }, numeric(nrow(tData)))
        
        
        pNAvalues[, vs] <- rowMins(pvals) * (NumReps + 1)
        qNAvalues[, vs] <- p.adjust(pNAvalues[, vs], method = "BH")
        setTxtProgressBar(pb, vs)
    }
    close(pb)
    pNAvalues[pNAvalues > 1] <- 1
    
    return(list(pNAvalues = pNAvalues, qNAvalues = qNAvalues))
}


#' Calculates the p-values and q-values for missing values in a data frame with
#' columns as samples and rows as features.
#'
#' @param Data A matrix containing the expression data with rows as features and
#' columns as samples
#' @param NumCond The number of conditions in the dataset
#' @param NumReps The number of replicates for each condition (needs to be same
#' for all conditions)
#'
#' @return A list containing the calculated p-values and q-values
#' (Benjamini-Hochberg) for missingness statistics
#'
#' @examples
#' Data <- matrix(rnorm(120), nrow = 10)
#' # Introduce some missingness
#' Data[sample(1:120, 40)] <- NA
#' NumCond <- 4
#' NumReps <- 3
#' res_misstest <- MissingStats(Data, NumCond, NumReps)
#' head(res_misstest$qNAvalues)
#' @importFrom matrixStats rowMins
#' @export
MissingStats <- function(Data, NumCond, NumReps) {
    Reps <- rep(seq_len(NumCond), NumReps)
    pNAvalues <- matrix(NA,
                        ncol = NumCond - 1, nrow = nrow(Data),
                        dimnames = list(
                            rows = rownames(Data),
                            cols = seq_len((NumCond - 1))
                        )
    )
    qNAvalues <- matrix(NA,
                        ncol = NumCond - 1, nrow = nrow(Data),
                        dimnames = list(
                            rows = rownames(Data),
                            cols = seq_len((NumCond - 1))
                        )
    )
    for (vs in 2:NumCond) {
        tData <- Data[, Reps == vs, drop=FALSE]
        trefData <- Data[, Reps == 1, drop=FALSE]
        tCompDat <- cbind(tData, trefData)
        qs <- quantile(tCompDat, probs = seq(0, 1, 0.01), na.rm = TRUE)
        
        pvals <- vapply(qs, function(q) {
            tCompDatQ <- tCompDat
            tCompDatQ[tCompDatQ < q] <- NA
            NAPDistr <- MissValPDistr(NumReps, sum(is.na(tCompDatQ)) /
                                          (nrow(tCompDatQ) * 2 * NumReps))
            statis <- rowSums(!is.na(tCompDatQ[, seq_len(NumReps)])) -
                rowSums(!is.na(tCompDatQ[, (NumReps + 1):(2 * NumReps)]))
            pvals_q <- NAPDistr[abs(statis) + 1]
            return(pvals_q)
        }, numeric(nrow(tCompDat)))
        pNAvalues[, vs - 1] <- rowMins(pvals) * (NumReps + 1)
        qNAvalues[, vs - 1] <- p.adjust(pNAvalues[, vs - 1], method = "BH")
    }
    
    pNAvalues[pNAvalues > 1] <- 1
    
    return(list(pNAvalues = pNAvalues, qNAvalues = qNAvalues))
}

#' Function to determine "optimal" fold-change and q-value thresholds using
#' a higher criticism method and an FC threshold based on a standard deviation
#' of one
#'
#' @param Pvalue A matrix of p-values for each condition
#' @param LogRatios A matrix of log ratios for each condition
#'
#' @return A vector containing the mean fold-change threshold and mean q-value
#' threshold
#'
#' @examples
#' Pvalue <- matrix(seq(0.01, 0.12, 0.01), nrow = 4, ncol = 3)
#' LogRatios <- matrix(c(
#'     1.2, 0.8, 1.5, -0.5, 0.2, 0.9, -1.1, 0.7, 1.8,
#'     -0.9, 0.3, 1.1
#' ), nrow = 4, ncol = 3)
#' thresholds <- FindFCandQlimAlternative(Pvalue, LogRatios)
#' print(thresholds)
#'
#' @export
#' @importFrom fdrtool hc.thresh
FindFCandQlimAlternative <- function(Pvalue, LogRatios) {
    BestComb <- c(0, 0)
    NumCond <- ncol(LogRatios)
    Pvalue[is.na(Pvalue)] <- 1
    
    NumTests <- ncol(Pvalue) / (NumCond)
    if (NumTests %% 1 != 0) {
        stop("Number of tests is not a multiple of the number of conditions")
    }
    if (NumTests != ncol(Pvalue) / (NumCond)) {
        stop("Number of tests is not a multiple of the number of conditions")
    }
    
    
    # Apply hc.thresh function to p-values for each column
    BestHCs <- apply(Pvalue, 2, function(pvals) {
        hc.thresh(pvals[pvals < 1], plot = FALSE)
    })
    
    # Calculate the standard deviation of LogRatios for each column
    BestFCs <- apply(LogRatios, 2, sd, na.rm = TRUE)
    
    BestHCs
    BestFCs
    
    # Calculate mean of all estimated thresholds
    
    return(c(mean(BestFCs), mean(BestHCs)))
}



#' Find Optimal Fold-Change and Q-Value Thresholds
#'
#' Identifies the optimal fold-change (FC) and q-value thresholds that maximize
#' the number of significant features identified across different statistical
#' tests. It processes a matrix of q-values, applying various fold-change
#' thresholds,
#' and computes the distribution of significant features for each q-value
#' threshold
#' to determine the optimal combination of thresholds.
#'
#' @param Qvalue A matrix of q-values obtained from statistical tests such as
#' PolySTest, limma, rank products, and permutation tests. NA values in the
#' matrix are treated as non-significant. The matrix should have the same number
#' of rows as the LogRatios matrix, and the number of columns should be a
#' multiple of the number of conditions, given by the number of different
#' statistical tests.
#' @param LogRatios A matrix of log2 fold-change ratios for the same
#' comparisons.
#'
#' @return A numeric vector containing two elements: the optimal fold-change
#' threshold and the optimal q-value threshold, which together maximize the
#' number of significant features detected.
#'
#' @details The function first replaces NA values in the Qvalue matrix with 1 to
#' denote non-significant results. It then identifies the smallest q-value
#' and calculates a range of q-values around this minimum. For each fold-change
#' threshold, the function adjusts the Qvalue matrix and calculates the number
#' of significant features for each q-value threshold. The combination yielding
#' the highest mean distribution of significant features is considered optimal.
#'
#' @examples
#' # Example Qvalue and LogRatios matrices
#' Qvalue <- matrix(seq(0.01, 0.12, 0.01), nrow = 4, ncol = 3)
#' LogRatios <- matrix(c(
#'     1.2, 0.8, 1.5, -0.5, 0.2, 0.9, -1.1, 0.7,
#'     1.8, -0.9, 0.3, 1.1
#' ), nrow = 4, ncol = 3)
#' # Find optimal thresholds
#' thresholds <- FindFCandQlim(Qvalue, LogRatios)
#' print(thresholds)
#'
#' @export
#' @importFrom matrixStats colMins
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @keywords internal
FindFCandQlim <- function(Qvalue, LogRatios) {
    message("Estimating suitable values for FDR and FC cutoff")
    BestComb <- c(0, 0)
    BestRegs <- 0
    NumCond <- ncol(LogRatios)
    
    Qvalue[is.na(Qvalue)] <- 1
    
    smallestq <- signif(min(Qvalue, na.rm = TRUE))
    qrange <- c(0.1, 0.2, 0.5) * 10^(rep(-10:0, each = 3))
    qrange <- qrange[which.min(abs(smallestq - qrange)):(length(qrange) - 2)]
    
    fcRange <- seq(0, max(abs(range(LogRatios, na.rm = TRUE))), length = 100)
    NumTests <- ncol(Qvalue) / (NumCond)
    if (NumTests %% 1 != 0) {
        stop("Number of tests is not a multiple of the number of conditions")
    }
    if (NumTests != ncol(Qvalue) / (NumCond)) {
        stop("Number of tests is not a multiple of the number of conditions")
    }
    
    NumThreads <- get_numthreads()
    cl <- parallel::makeCluster(NumThreads)
    parallel::clusterExport(cl,
                            varlist = c(
                                "Qvalue", "NumCond", "LogRatios",
                                "qrange", "NumTests"
                            ),
                            envir = environment()
    )
    
    
    BestVals <- parallel::parLapply(cl, fcRange, function(fc) {
        # Initialize variables to track the best regulatory value and combination
        localBestRegs <- 0
        localBestComb <- c(0, 0)
        
        for (t in seq_len(NumTests)) {
            # Modify Qvalue based on FC for current test
            tvals <- Qvalue[, (t - 1) * NumCond + seq_len(NumCond), drop = FALSE]
            tvals[LogRatios < fc & LogRatios > -fc] <- 1
            # Update Qvalue for current test
            Qvalue[, (t - 1) * NumCond + seq_len(NumCond)] <- tvals
            
            # Iterate over q-value range
            for (qlim in qrange) {
                # Calculate distribution of significant features for current
                # q-value threshold
                alldistr <- vapply(seq_len(NumCond), function(cond) {
                    distr <- table(rowSums(Qvalue[, (t - 1) * NumCond + cond,
                                                  drop = FALSE
                    ] < qlim, na.rm = TRUE))
                    sum(distr[-1], na.rm = TRUE) / sum(distr, na.rm = TRUE)
                }, numeric(1))
                
                if (mean(alldistr, na.rm = TRUE) > localBestRegs) {
                    localBestRegs <- mean(alldistr, na.rm = TRUE)
                    localBestComb <- c(fc, qlim)
                }
            }
            
        }
        c(localBestRegs, localBestComb)
    })
    parallel::stopCluster(cl)
    
    # Find the global best combination across all fold-change thresholds
    # Convert list to matrix for easier handling
    globalBest <- do.call(rbind, BestVals)
    bestRow <- which.max(globalBest[, 1])
    BestComb <- globalBest[bestRow, -1]
    BestRegs <- globalBest[bestRow, 1]
    
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
    message("Calculating PolySTest FDRs ...")
    UnifiedQvalue <- matrix(NA, ncol = NumComps, nrow = nrow(Qvalue))
    colnames(UnifiedQvalue) <- paste(
        "q_values_polystest_fdr_",
        seq_len(NumComps)
    )
    rownames(UnifiedQvalue) <- rownames(Qvalue)
    for (i in seq_len(NumComps)) {
        UnifiedQvalue[, i] <-
            colMins(apply(
                Qvalue[, seq(i, ncol(Qvalue), NumComps)], 1,
                p.adjust, "hommel"
            ), na.rm = TRUE)
    }
    UnifiedQvalue
}


########### Functions for CLI ############

## Parameter checks
validate_parameters <- function(pars) {
    validate_numeric_param(pars$numreps, "Number of replicates",
                           larger_than = 0, integer = TRUE
    )
    validate_numeric_param(pars$numcond, "Number of conditions",
                           larger_than = 0, integer = TRUE
    )
    validate_numeric_param(pars$refcond, "Reference condition",
                           larger_than = -1, integer = TRUE
    )
    validate_numeric_param(pars$firstquantcol, "First quant column",
                           larger_than = 0, integer = TRUE
    )
    validate_numeric_param(pars$threads, "Number of threads",
                           larger_than = 0, integer = TRUE
    )
    
    validate_choice_param(
        pars$normalization,
        "Normalization",
        c(
            "none", "quantile", "median", "mean",
            "sum_abs", "cyclicloess"
        )
    )
    validate_choice_param(pars$delim, "Delimiter", c("tab", ";", ",", "|"))
    validate_choice_param(pars$decimal, "Decimal", c(".", ","))
    validate_logical_param(pars$paired, "Paired")
    validate_logical_param(pars$header, "Header")
    validate_logical_param(pars$rep_grouped, "Replicate groups")
}

validate_numeric_param <- function(param, name, larger_than = 0,
                                   integer = FALSE) {
    if (integer && !is.integer(param)) {
        stop(sprintf("%s must be an integer.", name), call. = FALSE)
    }
    if (param <= larger_than) {
        stop(sprintf("%s must be a number > %d.", name, larger_than),
             call. = FALSE
        )
    }
}

validate_choice_param <- function(param, name, choices) {
    if (!(param %in% choices)) {
        warning(sprintf(
            "%s must be one of the following: %s. Defaulting to 'none'",
            name, paste(choices, collapse = ", ")
        ), call. = FALSE)
        return("none") # or return the first choice or any specific default
    }
    return(param)
}

validate_logical_param <- function(param, name) {
    if (!is.logical(param)) {
        stop(sprintf("%s parameter should be 'true' or 'false'", name),
             call. = FALSE
        )
    }
}


#' Update Conditions with Longest Common Subsequence
#'
#' This function iterates over a set of conditions and updates each condition
#' with the longest common subsequence of column names associated with that
#' condition.
#'
#' @param fulldata A SummarizedExperiment or derived object that contains
#' the quantitative data as specified for PolySTest
#' @param default A vector of the length of the number of conditions suggesting
#' their names
#'
#' @return SummarizedExperiment with updated conditions
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' 
#' @export
#' @examples
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays = list(count = matrix(rnorm(200),
#'     ncol = 10
#' )))
#' metadata(se) <- list(NumCond = 2, NumReps = 5)
#' rownames(colData(se)) <- paste0(
#'     rep(c("CondA_Rep", "CondB_Rep"), 5),
#'     rep(seq_len(5), each = 2)
#' )
#' default_conditions <- c("Condition_A", "Condition_B")
#' updated_conditions <- update_conditions_with_lcs(se, default_conditions)
#' print(colData(updated_conditions))
update_conditions_with_lcs <- function(fulldata, default = NULL) {
    # check for consistency
    check_for_polystest(fulldata)
    
    # extract data from SummarizedExperiement
    quant_cols <- colnames(SummarizedExperiment::assay(fulldata))
    numCond <- metadata(fulldata)$NumCond
    numReps <- metadata(fulldata)$NumReps
    if (is.null(default)) {
        default <- paste("C", seq_len(numCond), sep = "_")
    }
    conditions <- default
    
    # Function to find the longest common subsequence
    lcs <- function(a, b) {
        A <- strsplit(a, "")[[1]]
        B <- strsplit(b, "")[[1]]
        L <- matrix(0, length(A), length(B))
        ones <- which(outer(A, B, "=="), arr.ind = TRUE)
        ones <- ones[order(ones[, 1]), , drop = FALSE]
        if (nrow(ones) > 0) {
            for (i in seq_len(nrow(ones))) {
                v <- ones[i, , drop = FALSE]
                L[v] <- ifelse(any(v == 1), 1, L[v - 1] + 1)
            }
            return(paste0(
                A[(-max(L) + 1):0 + which(L == max(L),
                                          arr.ind = TRUE
                )[1]],
                collapse = ""
            ))
        } else {
            return(NA)
        }
    }
    
    
    # Iterate over each condition
    for (i in seq_len(numCond)) {
        condNames <- quant_cols[(0:(numReps - 1)) * numCond + i]
        lcsName <- condNames[1]
        if (length(condNames) > 1) {
            
            for (j in 2:length(condNames)) {
                if (!is.na(lcsName)) {
                    lcsName <- lcs(lcsName, condNames[j])
                }
            }
        }
        if (!is.na(lcsName) & !lcsName %in% conditions) {
            conditions[i] <- lcsName
        }
    }
    
    # Update colData of fulldata
    SummarizedExperiment::colData(fulldata)$Condition <- rep(conditions, numReps)
    
    return(fulldata)
}

#' Create All Pairwise Comparisons
#'
#' This function generates a matrix of all pairwise comparisons between the
#' provided conditions.
#' Optionally, a reference condition can be specified, and comparisons will be
#' made between the reference condition and all other conditions.
#'
#' @param conditions A character vector of condition names.
#' @param refCond An integer indicating the index of the reference condition
#' within the
#' `conditions` vector. If `refCond` is greater than 0, comparisons are made
#' between the reference condition and all other conditions. If `refCond` is 0,
#' all possible pairwise comparisons are made. Default is 0.
#'
#' @return A matrix with two columns, representing all possible pairwise
#' comparisons between the specified conditions. If a reference condition is
#' specified, it appears in the first column of every row.
#'
#' @examples
#' conditions <- c("Cond1", "Cond2", "Cond3")
#' # Generate all pairwise comparisons
#' allComps <- create_pairwise_comparisons(conditions, refCond = 0)
#' allComps
#'
#' # Generate comparisons with a specific condition as reference
#' refComps <- create_pairwise_comparisons(conditions, refCond = 1)
#' refComps
#'
#' @importFrom knitr kable
#' @export
create_pairwise_comparisons <- function(conditions, refCond) {
    if (refCond > 0) {
        allComps <- cbind(
            conditions[refCond],
            conditions[conditions != conditions[refCond]]
        )
    } else {
        allComps <- t(combn(conditions, 2))
    }
    colnames(allComps) <- c("Condition A", "Condition B")
    message("All pairwise comparison between conditions:")
    cat(knitr::kable(allComps), sep = "\n")
    return(allComps)
}

# Convert pairwise comparisons to numeric indices for statistical analysis

create_ratio_matrix <- function(fulldata, allComps) {
    # check for consistency
    check_for_polystest(fulldata)
    # extract data from SummarizedExperiement
    conditions <- unique(SummarizedExperiment::colData(fulldata)$Condition)
    NumCond <- metadata(fulldata)$NumCond
    NumReps <- metadata(fulldata)$NumReps
    dat <- SummarizedExperiment::assay(fulldata)
    
    # Make indices for pairings
    valComps <- matrix(NA, nrow = nrow(allComps), ncol = 2)
    valComps[, 1] <- match(allComps[, 1], conditions)
    valComps[, 2] <- match(allComps[, 2], conditions)    
    
    RR <- matrix(NA, ncol = nrow(allComps) * NumReps, nrow = 2)
    for (j in seq_len(nrow(allComps))) {
        RR[1, seq(j, nrow(allComps) * NumReps, nrow(allComps))] <-
            seq(valComps[j, 2], NumCond * NumReps, NumCond)
        RR[2, seq(j, nrow(allComps) * NumReps, nrow(allComps))] <-
            seq(valComps[j, 1], NumCond * NumReps, NumCond)
    }
    # Create ratio matrix
    
    MAData <- dat[, RR[1, ]] - dat[, RR[2, ]]
    
    rownames(MAData) <- rownames(dat)
    return(MAData)
}


# reducing repetition in limma test calls
fit_and_getvals <- function(lm.fitted) {
    
    lm.bayes <- limma::eBayes(lm.fitted)
    topTable(lm.bayes)
    plvalues <- lm.bayes$p.value
    qlvalues <- matrix(NA,
                       nrow = nrow(plvalues), ncol = ncol(plvalues),
                       dimnames = dimnames(plvalues)
    )
    plvalues[!is.finite(plvalues)] <- NA
    
    # qvalue correction
    for (i in seq_len(ncol(plvalues))) {
        p_vals <- na.omit(plvalues[, i])
        
        # Use tryCatch to attempt qvalue first and fall back to p.adjust with "BH" method
        tqs <- tryCatch({
            # Try qvalue method
            qvalue::qvalue(p_vals)$qvalues
        }, error = function(e) {
            message("qvalue failed, falling back to Benjamini-Hochberg: ", e$message)
            # Fallback to Benjamini-Hochberg method
            p.adjust(p_vals, method = "BH")
        })    
        
        qlvalues[names(tqs), i] <- tqs
    }
    
    return(list(
        plvalues = plvalues, qlvalues = qlvalues,
        Sds = sqrt(lm.bayes$s2.post)
    ))
}

# # Make this a hdden function as the arguments are too demanding
# ' Prepare Output Data for PolySTest Results
# '
# ' This function processes the results of PolySTest, including log-ratios,
# '  p-values, and q-values for each statistical test applied, and integrates
# '  them into the rowData of the provided SummarizedExperiment object. It
# '  optionally handles separate t-test q-values and unifies q-values across
# '  multiple tests.
# '
# ' @param fulldata A SummarizedExperiment object containing the initial dataset.
# ' @param Pvalue A matrix of p-values from the statistical tests with column
# ' names starting with "p_values_"
# ' @param Qvalue A matrix of q-values corresponding to the p-values with column
# ' names starting with "FDR_"
# ' @param LogRatios A matrix of log-ratio values for the comparisons with column
# ' names starting with "log_ratios_
# ' @param testNames A vector of names for each of the statistical tests
# ' performed.
# ' @param allComps A matrix specifying the pairs of conditions compared,
# '        each row represents a pair.
# '
# ' @details This function first checks for the presence of t-test results within
# '          the provided data, segregates them if present, and then optionally
# '          unifies q-values across tests if more than one test is specified
# '          (with Hommel correction for multiple testing).
# '          It then organizes and renames the matrices of log-ratios, p-values,
# '          and q-values according to the comparisons and tests performed, and
# '          merges these matrices into the rowData of the provided
# '          SummarizedExperiment object. Finally, it prints a summary of the
# '          number of features with FDR < 0.01 for each test and comparison.
# '
#' @return The updated SummarizedExperiment object with additional columns in
# '         rowData for log-ratios, p-values, and q-values.
# '
# ' @examples
# ' # Assuming 'fulldat' is your SummarizedExperiment object, 'Pvalue', 'Qvalue',
# ' # and 'LogRatios' are matrices of your test results, 'testNames' is your
# ' # vector of test names, and 'allComps' specifies your condition pairs:
# ' # fulldat <- prepare_output_data(fulldat, Pvalue, Qvalue, LogRatios,
# ' #                                testNames, allComps)
# ' data(liver_example)
# ' rdata <- rowData(liver_example)
# ' Pvalue <- rdata[, grep("p_value_", colnames(rdata))]
# ' Qvalue <- rdata[, grep("FDR_", colnames(rdata))]
# ' LogRatios <- rdata[, grep("log_ratios_", colnames(rdata))]
# ' prepare_output_data(liver_example, Pvalue, 
# '                     Qvalue, LogRatios,
# '                    c("perm_test", "limma"), 
# '                    c("FO.Rep.", "TTA.Rep."))
# '
# ' @importFrom SummarizedExperiment rowData
# ' 
# ' @importFrom knitr kable
# ' @importFrom S4Vectors metadata
# ' @export
prepare_output_data <- function(fulldata, Pvalue, Qvalue, LogRatios,
                                testNames, allComps) {
    num_tests <- length(testNames)
    numComps <- nrow(allComps)
    testNames2 <- testNames
    num_tests2 <- num_tests
    
    # Separate t-test p-values
    if ("t_test" %in% testNames) {
        ttestQvalue <- Qvalue[, grep("q_values_t_test", colnames(Qvalue))]
        Qvalue <- Qvalue[, -grep("q_values_t_test", colnames(Qvalue))]
        num_tests <- num_tests - 1
    }
    
    # Unify q-values if needed (assuming UnifyQvals is a function to unify
    # q-values)
    if (num_tests > 1) {
        message("Unifying q-values across tests ...")
        Qvalue <- cbind(UnifyQvals(Qvalue, numComps, num_tests), Qvalue)
        num_tests2 <- num_tests + 1
        testNames2 <- c("PolySTest", testNames)
    }
    
    if ("t_test" %in% testNames) {
        Qvalue <- cbind(Qvalue, ttestQvalue)
        # move "t_test" to end
        testNames2 <- c(testNames2[testNames2 != "t_test"], "t_test")
        num_tests <- num_tests + 1
    }
    
    # Set column names for log-ratios, p-values, and q-values
    compNames <- apply(allComps, 1, function(x) {
        paste(x[2], "vs", x[1],
              sep = "_"
        )
    })
    colnames(LogRatios) <- paste("log_ratios", compNames, sep = "_")
    colnames(Pvalue) <- paste("p_values", rep(testNames, each = numComps),
                              rep(compNames, num_tests),
                              sep = "_"
    )
    colnames(Qvalue) <- paste("FDR", rep(testNames2, each = numComps),
                              rep(compNames, num_tests2),
                              sep = "_"
    )
    
    # Combine all data into a single data frame
    SummarizedExperiment::rowData(fulldata) <- cbind(SummarizedExperiment::rowData(fulldata), LogRatios, Qvalue, Pvalue)
    
    # Assuming FullReg contains q-values and you are interested in features with
    # FDR < 0.01
    message("------- Summary of Results --------")
    message("Number of differentially regulated features with FDR < 0.01:")
    
    # Calculate the number of features with FDR < 0.01
    significantFeatures <- apply(Qvalue, 2, function(x) {
        sum(x < 0.01,
            na.rm = TRUE
        )
    })
    
    # Print the summary
    cat(
        knitr::kable(matrix(significantFeatures,
                            nrow = length(compNames),
                            dimnames = list(rows = compNames, cols = testNames2)
        )),
        sep = "\n"
    )
    
    # Adding test details to metadata
    S4Vectors::metadata(fulldata)$testNames <- testNames2
    S4Vectors::metadata(fulldata)$allComps <- allComps
    S4Vectors::metadata(fulldata)$compNames <- paste0(allComps[, 2], "_vs_", allComps[, 1])
    
    return(fulldata)
}

#' Check SummarizedExperiment for PolySTest Requirements
#'
#' Performs checks on a SummarizedExperiment object to ensure it is properly
#' formatted for PolySTest analysis. It checks for specific metadata properties,
#' the number of assays, and the distribution of conditions and replicates.
#'
#' @param se A SummarizedExperiment object.
#'
#' @return Invisible TRUE if checks pass; otherwise, warnings or errors are
#' thrown.
#'
#' @examples
#' data(liver_example)
#' check_for_polystest(liver_example)
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' 
#' @export
check_for_polystest <- function(se) {
    # Check for required metadata
    if (!("NumReps" %in% names(metadata(se))) ||
        !("NumCond" %in% names(metadata(se)))) {
        stop("Metadata must contain 'NumReps' (number of replicates) and
        'NumCond' (number of experimental conditions).")
    }
    
    # Check values of NumReps and NumCond
    if (metadata(se)$NumReps <= 1) {
        warning("NumReps is 1, which may not be suitable for analysis.")
    }
    if (metadata(se)$NumCond <= 1) {
        warning("NumCond is 1, which may not be suitable for analysis.")
    }
    
    # Check for only one assay
    if (length(SummarizedExperiment::assays(se)) != 1) {
        stop("SummarizedExperiment should contain only one assay.")
    }
    
    # Check if the number of columns in the assay is NumReps * NumCond
    expectedCols <- metadata(se)$NumReps * metadata(se)$NumCond
    if (ncol(assay(se)) != expectedCols) {
        stop("The number of columns in the assay does not match NumReps *
             NumCond.")
    }
    
    # Check if each condition has the same number of replicates
    condReps <- table(SummarizedExperiment::colData(se)$Condition)
    if (any(condReps != metadata(se)$NumReps)) {
        warning("Not all conditions have the same number of replicates.")
    }
    
    # Any other checks can be added here
    
    invisible(TRUE)
}


filterFC <- function(rdat, NumTests, NumComps, fclim = c(0, 0)) {
    Qvalue <- as.data.frame(rdat[, grep("^FDR", colnames(rdat))])
    if (ncol(Qvalue) > 0) {
        LogRatios <- as.data.frame(rdat[, grep("^log_ratios", colnames(rdat))])
        
        FCRegs <- Qvalue
        
        for (t in seq_len(NumTests)) {
            tsign <- Qvalue[, (t - 1) * (NumComps) + (seq_len(NumComps)),
                            drop = FALSE
            ]
            tsign[is.na(tsign)] <- 1
            tsign[LogRatios > fclim[1] & LogRatios < fclim[2]] <- 1
            FCRegs[, (t - 1) * (NumComps) + (seq_len(NumComps))] <- tsign
        }
        FCRegs
    } else {
        NULL
    }
}

#' Check Statistical Test and Comparison Names
#'
#' Verifies if the provided test names and comparison names are correct
#' and have been executed within the given `SummarizedExperiment` object.
#' It checks for the presence of specific metadata related to the tests
#' and comparisons to ensure that the requested analyses have been carried out.
#'
#' @param fulldata A `SummarizedExperiment` object containing the dataset
#' and metadata for statistical analyses.
#' @param compNames A character vector of comparison names to be verified
#' against the `SummarizedExperiment` metadata. If set to "all",
#' the function checks for the presence of any comparison names in the metadata.
#' @param testNames A character vector of statistical test names to be verified
#' against the `SummarizedExperiment` metadata.
#'
#' @details This function is used to ensure that the requested statistical tests
#' and comparisons have been correctly named and executed prior to further
#' analysis. It verifies that the names provided match those stored within the
#' metadata of the `SummarizedExperiment` object. If the names do not match or
#' the necessary metadata is missing, the function stops execution and returns
#' an error message indicating the issue.
#'
#' @return The function return the updated comparison names if the checks pass.
#'
#' @examples
#' data(liver_example)
#' compNames <- "all"
#' testNames <- c("limma", "t_test")
#' check_stat_names(liver_example, compNames, testNames)
#' 
#' @export
check_stat_names <- function(fulldata, compNames, testNames) {
    # Check for correct test and column names
    if (is.null(metadata(fulldata)$testNames) ||
        !all(testNames %in% metadata(fulldata)$testNames)) {
        stop("The test names are not correct or the tests have not been carried
         out.")
    }
    
    # Check for correct comparison names
    if (is.null(metadata(fulldata)$compNames)) {
        stop("The comparison names are not correct or the comparisons have not
        been carried out.")
    }
    if (compNames[1] == "all") {
        compNames <- metadata(fulldata)$compNames
    }
    if (!all(compNames %in% metadata(fulldata)$compNames)) {
        stop("The comparison names are not correct or the comparisons have not
        been carried out.")
    }
    
    return(compNames)
}

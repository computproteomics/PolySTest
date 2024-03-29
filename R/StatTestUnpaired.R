#' PolySTest for unpaired tests
#'
#' Combining the power of different statistical tests
#'
#' @param fulldata A SummarizedExperiment or derived object that contains
#' the quantitative data as required for PolySTest
#' @param allComps A matrix containing the reference matrix specifying the
#' pairs of conditions to compare (each comparison given as separate row)
#' @param statTests A character vector specifying the statistical tests to be
#' used. The available tests are: "limma", "Miss_Test", "t-test",
#' "rank_products", and "permutation_test"
#'
#' @details This function performs unpaired statistical tests on the data in
#' 'fulldata' using the pairs of conditions specified in 'allComps'. It
#' calculates the p-values and q-values for each row for the statistical tests
#' used. The statistical tests available are: limma, Miss_Test, t-test,
#' rank_products, and a permutation_test based on t values. The function returns
#' a SummarizedExperiment object with added columns for p-values and q-values in
#' rowData.
#'
#' @return SummarizedExperiment with added columns for p-values and q-values
#' in rowData
#'
#' @examples
#' # Creating mock quantitative data and sample metadata
#' library(SummarizedExperiment)
#' quantData <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#' colnames(quantData) <- c(
#'     paste("Sample", seq_len(5), "_Condition_A", sep = ""),
#'     paste("Sample", seq_len(5), "_Condition_B", sep = "")
#' )
#' rownames(quantData) <- paste("Gene", seq_len(200))
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
#' # Running the PolySTest_unpaired function
#' results <- PolySTest_unpaired(fulldata, allComps)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#'
PolySTest_unpaired <- function(fulldata, allComps,
                               statTests = c(
                                   "limma", "Miss_Test", "t_test",
                                   "rank_products",
                                   "permutation_test"
                               )) {
    # check fulldata
    check_for_polystest(fulldata)

    # Extracting the assay data
    Data <- assay(fulldata, "quant")

    # Extracting the metadata
    NumReps <- metadata(fulldata)$NumReps
    NumCond <- metadata(fulldata)$NumCond

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

    # Normalize row-wise by mean
    Data <- Data - rowMeans(Data, na.rm = TRUE)

    # Prepare output data
    tests <- statTests
    # Check for the right test names
    if (any(!tests %in% c(
        "limma", "Miss_Test", "t_test",
        "rank_products", "permutation_test"
    ))) {
        stop("Invalid test name(s) specified. They should one or more of:
         'limma', 'Miss_Test', 't_test', 'rank_products', 'permutation_test'")
    }
    p_values <- q_values <- matrix(NA,
        nrow = nrow(Data),
        ncol = length(tests) * NumComps
    )
    rownames(p_values) <- rownames(q_values) <- rownames(Data)
    colnames(p_values) <- paste0(
        "p_values_", rep(tests, each = NumComps),
        "_", rep(seq_len(NumComps), length(tests))
    )
    colnames(q_values) <- paste0(
        "q_values_", rep(tests, each = NumComps),
        "_", rep(seq_len(NumComps), length(tests))
    )

    ## limma
    Sds <- NULL
    if (any("limma" %in% tests)) {
        message("Running limma tests")
        lm_out <- limma_unpaired(Data, NumCond, NumReps, RRCateg)
        p_values[, grep("p_values_limma", colnames(p_values))] <-
            lm_out$plvalues
        q_values[, grep("q_values_limma", colnames(q_values))] <-
            lm_out$qlvalues
        Sds <- lm_out$Sds
        message("limma completed")
    }

    if (any("Miss_Test" %in% tests)) {
        message("Running Miss test")
        MissingStats <- MissingStatsDesign(Data, RRCateg, NumCond, NumReps)
        p_values[, grep("p_values_Miss_Test", colnames(p_values))] <-
            MissingStats$pNAvalues
        q_values[, grep("q_values_Miss_Test", colnames(q_values))] <-
            MissingStats$qNAvalues
        message("Miss test completed")
    }

    ## rank products + t-test
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
    pb <- txtProgressBar(0.9, NumComps)

    for (vs in seq_len(NumComps)) {
        if (!is.null(shiny::getDefaultReactiveDomain())) {
            shiny::setProgress(0.1 + 0.3 / (NumComps) * vs,
                detail = paste(
                    "tests for comparison", vs,
                    "of", NumComps
                )
            )
        }
        tData <- Data[, Reps == RRCateg[1, vs]]
        trefData <- Data[, Reps == RRCateg[2, vs]]

        ## t-test
        if (any("t_test" %in% tests)) {
            ttest_out <- ttest_unpaired(tData, trefData)
            p_values[, grep("p_values_t_test", colnames(p_values))[vs]] <-
                ttest_out$ptvalues
            q_values[, grep("q_values_t_test", colnames(q_values))[vs]] <-
                ttest_out$qtvalues
        }

        ## rank products
        if (any("rank_products" %in% tests)) {
            rp_out <- rp_unpaired(tData, trefData)
            p_values[, grep(
                "p_values_rank_products",
                colnames(p_values)
            )[vs]] <-
                rp_out$pRPvalues
            q_values[, grep(
                "q_values_rank_products",
                colnames(q_values)
            )[vs]] <-
                rp_out$qRPvalues
        }

        ## Permutation tests
        if (any("permutation_test" %in% tests)) {
            perm_out <- perm_unpaired(tData, trefData)
            p_values[, grep(
                "p_values_permutation_test",
                colnames(p_values)
            )[vs]] <-
                perm_out$pPermutvalues
            q_values[, grep(
                "q_values_permutation_test",
                colnames(q_values)
            )[vs]] <-
                perm_out$qPermutvalues
        }

        lratios <- cbind(
            lratios,
            rowMeans(Data[, Reps == RRCateg[1, vs]], na.rm = TRUE) -
                rowMeans(Data[, Reps == RRCateg[2, vs]], na.rm = TRUE)
        )
        setTxtProgressBar(pb, vs)
    }
    close(pb)

    # Prepare output data
    fulldata <- prepare_output_data(
        fulldata, p_values, q_values, lratios,
        tests, allComps
    )

    return(fulldata)
}

#' Perform unpaired limma analysis
#'
#' This function performs unpaired limma analysis on Data.
#' @param Data A matrix of quantitative expression data.
#' @param NumCond The number of conditions in the experiment.
#' @param NumReps The number of replicates per condition.
#' @param RRCateg A matrix specifying the conditons to be compared
#' (each comparison given as separate column).
#' @return A list containing the following results:
#'  - plvalues: The p-values from limma tests.
#'  - qlvalues: The q-values from limma tests.
#'  - Sds: The standard deviations of the Bayesian linear model.
#'  @details This function performs unpaired limma analysis on Data. It
#'  calculates the p-values and q-values for each row, indicating the
#'  significance of the difference between the two datasets.RRCateg is a
#'  matrix specifying the conditons to be compared (each comparison given
#'  as separate column). It thus has always 2 rows and n columns, where n
#'  is the number of comparisons. The function returns a list containing
#'  the p-values and q-values for each comparison, as well as the standard
#'  deviations of the Bayesian linear model.
#' @keywords limma unpaired analysis
#' @export
#' @import limma
#' @import qvalue
#' @examples
#' dataMatrix <- matrix(rnorm(900), ncol = 9)
#' NumCond <- 3
#' NumReps <- 3
#' colnames(dataMatrix) <- rep(c("A", "B", "C"), each = 3)
#' #  Specifying comparisons
#' RRCateg <- matrix(c(1, 2, 2, 3), nrow = 2, ncol = 2)
#' # Run function
#' results <- limma_unpaired(dataMatrix, NumCond, NumReps, RRCateg)
#' print(results$plvalues)
#'
limma_unpaired <- function(Data, NumCond, NumReps, RRCateg) {
    Reps <- rep(seq_len(NumCond), NumReps)
    NumComps <- ncol(RRCateg)
    design <- model.matrix(~ 0 + factor(Reps - 1))
    colnames(design) <- paste("i", c(seq_len(NumCond)), sep = "")
    contrasts <- NULL
    First <- 1
    for (i in seq_len(NumComps)) {
        contrasts <- append(contrasts, paste(colnames(design)[RRCateg[2, i]],
            "-", colnames(design)[RRCateg[1, i]],
            sep = ""
        ))
    }
    contrast.matrix <- limma::makeContrasts(
        contrasts = contrasts,
        levels = design
    )
    lm.fitted <- limma::lmFit(Data, design)

    lm.contr <- limma::contrasts.fit(lm.fitted, contrast.matrix)
    lm.bayes <- limma::eBayes(lm.contr)
    topTable(lm.bayes)
    plvalues <- lm.bayes$p.value
    qlvalues <- matrix(NA,
        nrow = nrow(plvalues), ncol = ncol(plvalues),
        dimnames = dimnames(plvalues)
    )
    # qvalue correction
    for (i in seq_len(ncol(plvalues))) {
        tqs <- qvalue::qvalue(na.omit(plvalues[, i]))$qvalues
        qlvalues[names(tqs), i] <- tqs
    }

    return(list(
        plvalues = plvalues, qlvalues = qlvalues,
        Sds = sqrt(lm.bayes$s2.post)
    ))
}

#' Perform unpaired t-tests on two datasets
#'
#' @details
#' This function performs unpaired t-tests between corresponding rows of two
#' datasets.
#' It calculates the p-values and q-values for each row, indicating the
#' significance of the difference between the two datasets. We require
#'  providing the same number of samples (columns) per group.
#'
#' @param tData A matrix or data frame with the quantitative features (via rows)
#' of the first group
#' @param trefData A matrix or data frame with the quantitative features
#' (via rows) of the second group
#'
#' @return A list containing the p-values and q-values for each row
#'
#' @examples
#' tData <- matrix(rnorm(1000), nrow = 100)
#' trefData <- matrix(rnorm(1000), nrow = 100)
#' result <- ttest_unpaired(tData, trefData)
#' print(result$ptvalues)
#' print(result$qtvalues)
#'
#' @export
#' @import qvalue
ttest_unpaired <- function(tData, trefData) {
    # Check for the same column numbers
    if (ncol(tData) != ncol(trefData)) {
        stop("The number of columns in the two datasets must be the same")
    }

    ## t-tests
    tptvalues <- vapply(seq_len(nrow(tData)), function(pep) {
        ifelse(sum(!is.na(tData[pep, ])) > 1 & sum(!is.na(trefData[pep, ])) > 1,
            t.test(unlist(tData[pep, ]), unlist(trefData[pep, ]))$p.value,
            NA_real_
        )
    }, numeric(1))
    names(tptvalues) <- rownames(tData)
    ptvalues <- tptvalues
    tqs <- qvalue::qvalue(na.omit(ptvalues))$qvalues
    qtvalues <- rep(NA, length(ptvalues))
    names(qtvalues) <- names(ptvalues)
    qtvalues[names(tqs)] <- tqs

    return(list(ptvalues = ptvalues, qtvalues = qtvalues))
}


#' Perform unpaired rank products test
#'
#' @details
#' This function calculates the p-values and q-values using rankd product
#' statistics. The function requires having the same number of samples per
#' group. The function uses parallel computing to speed up the calculations.
#'
#' @param tData The data matrix for the test group (features are rows).
#' @param trefData The data matrix for the reference group (features are rows).
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
    # calculate NumRPPairs random pairing combinations and then take mean of
    # p-values
    if (exists("NumRPPairs")) {
        NumRPPairs <- NumRPPairs
    } else {
        NumRPPairs <- 100
    }

    # Check for the same column numbers
    if (ncol(tData) != ncol(trefData)) {
        stop("The number of columns in the two datasets must be the same")
    }


    NumReps <- ncol(tData)
    tpRPvalues <- matrix(NA,
        ncol = NumRPPairs, nrow = nrow(tData),
        dimnames = list(rows = rownames(tData), cols = seq_len(NumRPPairs))
    )
    NumThreads <- get_numthreads()
    cl <- parallel::makeCluster(NumThreads)
    parallel::clusterExport(
        cl = cl,
        varlist = c(
            "NumReps", "tData", "trefData",
            "RPStats"
        ), envir = environment()
    )
    #parallel::clusterEvalQ(cl = cl)

    RPparOut <- parallel::parLapply(cl, seq_len(NumRPPairs), function(x) {
        tRPMAData <- tData[, sample(seq_len(NumReps))] -
            trefData[, sample(seq_len(NumReps))]
        # Up
        RPMAUp_pvalues <- RPStats(tRPMAData, NumReps)
        # Down
        RPMADown_pvalues <- RPStats(-tRPMAData, NumReps)
        ttt <- rowMins(cbind(RPMAUp_pvalues, RPMADown_pvalues), na.rm = TRUE) * 2
        ttt[ttt > 1] <- 1
        names(ttt) <- names(RPMAUp_pvalues)
        ttt
    })
    stopCluster(cl)

    for (p in seq_len(NumRPPairs)) {
        # print(RPparOut[[p]])
        tpRPvalues[names(RPparOut[[p]]), p] <- RPparOut[[p]]
    }
    tpRPvalues[!is.finite(tpRPvalues)] <- NA
    pRPvalues <- rowMeans(tpRPvalues, na.rm = TRUE)
    qRPvalues <- rep(NA, length(pRPvalues))
    names(qRPvalues) <- names(pRPvalues)
    tqs <- p.adjust(na.omit(pRPvalues), method = "BH")
    qRPvalues[names(tqs)] <- tqs

    return(list(pRPvalues = pRPvalues, qRPvalues = qRPvalues))
}

#' Perform unpaired permutation tests
#'
#' This function performs permutation testing for unpaired data. The permutation
#' testing is based on comparing the t-values of the real data with the t-values
#' of the permuted data.
#'
#' @param tData The data matrix for the test group (features are rows).
#' @param trefData The data matrix for the reference group (features are rows).
#'
#' @return A list containing the p-values and q-values for the permutation test.
#'
#' @details The function adds columns from the randomized full set to reach the
#' minimum number of permutation columns (NumPermCols) replicates. It
#' randomizes the sign as well to avoid tendencies to one or the other side.
#' In the unpaired case, it also normalizes by the mean of the entire sample to
#' avoid strange effects. The function then performs permutation
#' testing using parallel computing, and calculates the p-values and q-values
#' based on the permutation results.
#' Both groups needs to consist of the same number of samples (columns).
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
    # Check for the same column numbers
    if (ncol(tData) != ncol(trefData)) {
        stop("The number of columns in the two datasets must be the same")
    }

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
    # In the unpaired case, also normalize by mean of the entire sample to avoid
    # strange effects
    tData <- tData - mean(as.numeric(unlist(tData)), na.rm = TRUE)
    trefData <- trefData - mean(as.numeric(unlist(trefData)), na.rm = TRUE)
    if (ncol(tData) * 2 < NumPermCols) {
        AddDat <- matrix(
            sample(as.vector(unlist(tData)),
                (NumPermCols - ncol(tData)) * nrow(tData),
                replace = TRUE
            ),
            nrow = nrow(tData)
        )
        PermData <- cbind(tData, AddDat)
        AddDat <- matrix(
            sample(as.vector(unlist(trefData)),
                (NumPermCols - ncol(trefData)) *
                    nrow(trefData),
                replace = TRUE
            ),
            nrow = nrow(trefData)
        )
        PermFullData <- cbind(PermData, trefData, AddDat)
    } else {
        PermFullData <- cbind(tData, trefData)
    }
    RealStats <- StatsForPermutTest(as.matrix(cbind(trefData, tData)),
        Paired = FALSE
    )
    # print(head(PermFullData))

    NumThreads <- get_numthreads()
    cl <- makeCluster(NumThreads)
    clusterExport(
        cl = cl, varlist = c(
            "NumReps", "PermFullData",
            "RPStats", "StatsForPermutTest"
        ),
        envir = environment()
    )
    #clusterEvalQ(cl = cl)
    PermutOut <- parallel::parLapply(cl, seq_len(NTests), function(x) {
        indat <- apply(PermFullData, 1, function(y) {
            sample(y, NumReps * 2) *
                sample(c(1, -1), NumReps * 2, replace = TRUE)
        })
        StatsForPermutTest(t(indat), FALSE)
    })
    stopCluster(cl)

    PermutOut <- matrix(unlist(PermutOut), nrow = nrow(tData))
    PermutOut[!is.finite(PermutOut)] <- NA
    RealStats[!is.finite(RealStats)] <- NA
    pPermutvalues <- apply(
        cbind(RealStats, PermutOut), 1,
        function(x) {
            ifelse(is.na(x[1]) | sum(!is.na(x)) == 0,
                NA,
                (1 + sum(x[1] < x[-1],
                    na.rm = TRUE
                )) /
                    (sum(!is.na(x)))
            )
        }
    )
    qPermutvalues <- rep(NA, length(pPermutvalues))
    names(qPermutvalues) <- names(pPermutvalues)
    tqs <- qvalue::qvalue(na.omit(pPermutvalues))$qvalues
    qPermutvalues[names(tqs)] <- tqs

    return(list(pPermutvalues = pPermutvalues, qPermutvalues = qPermutvalues))
}

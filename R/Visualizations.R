#' Set Graphics Layout
#'
#' Automatically sets the `mfrow` parameter for `par()` based on the total
#' number of plots and the maximum number of columns desired.
#'
#' @param num_total The total number of plots to display.
#' @param max_col The maximum number of columns for the layout.
#' 
#' @return Sets the `mfrow` parameter for `par()`.
#'
#' @examples
#' # This will set the layout to 2 rows of 3 columns
#' set_mfrow(num_total = 6, max_col = 3)
#' for (i in 1:6) hist(1:10)
#' par(mfrow = c(1, 1))
#' @export
#'
set_mfrow <- function(num_total, max_col) {
    if (num_total <= max_col) {
        par(mfrow = c(1, num_total))
    } else {
        par(mfrow = c(ceiling(num_total / max_col), max_col))
    }
}


#' Plot P-Value Distributions
#'
#' Generates histograms of p-value distributions for each test and comparison.
#'
#' @param fulldata A `SummarizedExperiment` object containing the dataset and
#' FDR/q-values from PolySTest.
#' @param compNames A character vector of comparison names. "all" selects all
#' comparisons.
#' @param testNames A character vector of test names used in the analysis.
#' Default values are "PolySTest", "limma", "Miss test", "rank products",
#' "permutation test", and "t-test".
#' @param testCols A character vector of colors for each test. Defaults to
#' c("#33AAAA", "#33AA33", "#AA3333", "#AA33AA", "#AAAA33", "#3333AA").
#' @param ... Additional arguments passed to `hist()`.
#' 
#' @return Creates histograms of p-value distributions for the specified tests
#'
#' @examples
#' # Assuming `fulldata` is a properly prepared `SummarizedExperiment` object
#' data(liver_example)
#' plotPvalueDistr(liver_example,
#'     compNames = c("HF.Rep._vs_TTA.Rep."),
#'     testCols = rainbow(5)
#' )
#' @export
plotPvalueDistr <- function(fulldata, compNames = "all",
                            testNames = c(
                                "limma", "Miss Test", "rank products",
                                "permutation test", "t-test"
                            ),
                            testCols = c(
                                "#33AAAA", "#33AA33", "#AA3333", "#AA33AA",
                                "#AAAA33", "#3333AA"
                            ), ...) {
    check_for_polystest(fulldata)

    # check compnames and testnames with metadata
    testNames2 <- make.names(testNames)
    testNames2 <- gsub("\\.", "_", testNames2)
    compNames <- check_stat_names(fulldata, compNames, testNames2)

    # Get the p-values from the right columns
    message("Plotting p-values")
    rdat <- SummarizedExperiment::rowData(fulldata)
    PValue <- as.matrix(rdat[, grep("^p_values_", colnames(rdat)),
        drop = FALSE
    ])
    # get only the p-values for the tests and comparisons
    PValue <- PValue[, grep(paste(paste0("^p_values_", testNames2),
        collapse = "|"
    ), colnames(PValue)),
    drop = FALSE
    ]
    PValue <- PValue[, grep(
        paste(paste0(compNames, "$"), collapse = "|"),
        colnames(PValue)
    ), drop = FALSE]

    NumTests <- length(testNames)
    NumComps <- length(compNames)
    set_mfrow((NumTests) * NumComps, NumTests)
    if (ncol(PValue) > 0) {
        for (i in seq_len(NumComps)) {
            for (j in seq_len(NumTests)) {
                hist(PValue[, (NumComps) * (j - 1) + i], 100,
                    main = testNames[j],
                    sub = compNames[i], col = testCols[j + 1],
                    xlab = "p-value", border = NA, ...
                )
            }
        }
        message("Plotting p-values finished")
    }
}

#' Plot Volcano Plots for PolySTest results
#'
#' This function creates volcano plots for all specified statistical tests and
#' comparisons using data from a `SummarizedExperiment` object. It highlights
#' selected proteins and applies fold-change and q-value limits for
#' visualization.
#'
#' @param fulldata A `SummarizedExperiment` object containing the data.
#' @param compNames A character vector of comparison names.
#' @param testNames A character vector of test names including "PolySTest",
#' "limma", "Miss test", "rank products", "permutation test", and "t-test".
#' @param sel_prots A numeric vector indicating selected features to be
#' visualized differently or "all" to select all features. Default is "all".
#' @param qlim A numeric value setting the q-value limit for the plots.
#' Default is 0.05.
#' @param fclim A numeric vector of length two setting the fold-change limits
#' for the plots. Default is c(0,0).
#' @param testCols A character vector of colors for each test. Default is a
#' predefined set of colors.
#' @param ... Additional arguments passed to the plot function.
#'
#' @return Creates volcano plots for the specified tests and comparisons.
#'
#' @examples
#' data(liver_example)
#' compNames <- c("HF.Rep._vs_TTA.Rep.")
#' plotVolcano(liver_example, compNames)
#'
#' @export
plotVolcano <- function(fulldata, compNames = "all",
                        testNames = c(
                            "PolySTest", "limma", "Miss Test",
                            "rank products", "permutation test", "t-test"
                        ),
                        sel_prots = "all", qlim = 0.05, fclim = c(0, 0),
                        testCols = c(
                            "#33AAAA", "#33AA33", "#AA3333", "#AA33AA",
                            "#AAAA33", "#3333AA"
                        ), ...) {
    check_for_polystest(fulldata)
    rdat <- SummarizedExperiment::rowData(fulldata)

    # check compnames and testnames with metadata
    testNames2 <- make.names(testNames)
    testNames2 <- gsub("\\.", "_", testNames2)
    compNames <- check_stat_names(fulldata, compNames, testNames2)

    Qvalue <- as.matrix(rdat[, grep("^FDR_", colnames(rdat)), drop = FALSE])
    LogRatios <- as.matrix(rdat[, grep("^log_ratios_", colnames(rdat)),
        drop = FALSE
    ])
    # get only the p-values for the tests and comparisons
    Qvalue <- Qvalue[, grep(
        paste(paste0("^FDR_", testNames2), collapse = "|"),
        colnames(Qvalue)
    ), drop = FALSE]
    LogRatios <- LogRatios[, grep(paste(paste0("^log_ratios_", compNames),
        collapse = "|"
    ), colnames(LogRatios)),
    drop = FALSE
    ]
    Qvalue <- Qvalue[, grep(
        paste(paste0(compNames, "$"), collapse = "|"),
        colnames(Qvalue)
    ), drop = FALSE]

    if (ncol(Qvalue) > 0) {
        message("Plotting volcano plots")
        NumTests <- length(testNames)
        NumComps <- length(compNames)
        set_mfrow(NumTests * NumComps, NumTests)
        if (all(sel_prots == "all") & length(sel_prots) > 0) {
            sel_prots <- seq_len(nrow(Qvalue))
        }
        for (i in seq_len(NumComps)) {
            for (j in seq_len(NumTests)) {
                plot(LogRatios[, i], -log10(Qvalue[, (NumComps) * (j - 1) + i]),
                    main = testNames[j], sub = compNames[i],
                    xlab = "log fold-change", ylab = "-log10(q)",
                    cex = colSelected(0.5, nrow(Qvalue), sel_prots, 1),
                    col = colSelected(
                        adjustcolor(testCols[j], alpha.f = 0.3), nrow(Qvalue),
                        sel_prots, "#FF9933"
                    ), pch = 16,
                    ylim = -log10(c(1, min(Qvalue, na.rm = TRUE))),
                    ...
                )

                abline(h = -log10(qlim), col = "#AA3333", lwd = 2)
                abline(v = fclim, col = "#AA3333", lwd = 2)
            }
        }
        message("Plotting volcano plots finished")
    }
}

#' Plot Expression Profiles
#'
#' This function plots expression profiles for selected features across
#' different conditions and comparisons. It supports both scaling and unscaled
#' profiles. It adds a circular plot to compare the different statistical tests
#'
#' @param fulldata A `SummarizedExperiment` object containing the data.
#' @param compNames A character vector of comparison names. "all" selects all
#' comparisons.
#' @param testNames A character vector of test names used in the analysis.
#' Default values are "PolySTest", "limma", "Miss test", "rank products",
#' "permutation test", and "t-test".
#' @param sel_prots A numeric vector with the indices of the selected features.
#'  Default is "all". These will still be filterd
#' @param profiles_scale Logical indicating if profiles should be scaled.
#' Default is TRUE.
#' @param qlim A numeric value indicating the q-value limit for significance.
#' @param fclim A numeric vector of length 2 indicating fold-change limits.
#'
#' @return Plots expression profiles for the selected features.
#' @examples
#' data(liver_example)
#' compNames <- c("HF.Rep._vs_TTA.Rep.")
#' plotExpression(liver_example)
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom gplots plotCI redblue
#' @import circlize
#'
plotExpression <- function(fulldata, compNames = "all",
                           testNames = c(
                               "PolySTest", "limma", "Miss Test",
                               "rank products", "permutation test", "t-test"
                           ),
                           sel_prots = "all", profiles_scale = TRUE,
                           qlim = 0.05, fclim = c(0, 0)) {
    check_for_polystest(fulldata)
    rdat <- SummarizedExperiment::rowData(fulldata)

    par(mfrow = c(1, 3))

    # check compnames and testnames with metadata
    testNames2 <- make.names(testNames)
    testNames2 <- gsub("\\.", "_", testNames2)
    compNames <- check_stat_names(fulldata, compNames, testNames2)

    NumComps <- length(compNames)
    NumTests <- length(testNames)

    Qvalue <- as.matrix(rdat[, grep("^FDR_", colnames(rdat)), drop = FALSE])
    LogRatios <- as.matrix(rdat[, grep("^log_ratios_", colnames(rdat))],
        drop = FALSE
    )
    # get only the p-values for the tests and comparisons
    Qvalue <- Qvalue[, grep(
        paste(paste0("^FDR_", testNames2), collapse = "|"),
        colnames(Qvalue)
    ), drop = FALSE]
    LogRatios <- LogRatios[, grep(paste(paste0("^log_ratios_", compNames),
        collapse = "|"
    ), colnames(LogRatios)),
    drop = FALSE
    ]
    Qvalue <- Qvalue[, grep(
        paste(paste0(compNames, "$"), collapse = "|"),
        colnames(Qvalue)
    ), drop = FALSE]

    dat <- SummarizedExperiment::assay(fulldata)
    rownames(Qvalue) <- rownames(LogRatios) <- rownames(rdat)
    if (all(sel_prots == "all")) {
        sel_prots <- seq_len(nrow(Qvalue))
    }

    if (ncol(Qvalue) > 0 & length(sel_prots) > 1) {
        message("plotting expression profiles")
        rdat <- rowData(fulldata)
        NumCond <- metadata(fulldata)$NumCond
        NumReps <- metadata(fulldata)$NumReps
        FCRegs <- filterFC(rdat, NumTests, NumComps, fclim)

        # CI plots of max 30 features
        SubSetQval <- Qvalue[sel_prots, , drop = FALSE]
        SubSetLR <- LogRatios[sel_prots, , drop = FALSE]
        SubSetLR <- SubSetLR[
            order(rowMins(
                SubSetQval[, seq_len(NumComps),
                    drop = FALSE
                ],
                na.rm = TRUE
            )), ,
            drop = FALSE
        ]

        SubSet <- SubSetLR[seq_len(min(nrow(SubSetLR), 30)), , drop = FALSE]
        indices <- rownames(SubSet)
        tdat <- as.matrix(dat[rownames(SubSet),
            (rep(seq_len(NumReps), NumCond) - 1) * NumCond +
                rep(seq_len(NumCond), each = NumReps),
            drop = FALSE
        ])
        rownames(tdat) <- strtrim(rownames(tdat), 20)
        MeanSet <- SDSet <- matrix(NA,
            nrow = nrow(tdat), ncol = NumCond, dimnames =
                list(x = rownames(tdat), y = paste(
                    "Condition",
                    seq_len(NumCond)
                ))
        )
        for (c in seq_len(NumCond)) {
            MeanSet[, c] <- rowMeans(tdat[, seq_len(NumReps) + (c - 1) *
                NumReps,
            drop = FALSE
            ], na.rm = TRUE)

            SDSet[, c] <- rowSds(tdat[, seq_len(NumReps) + (c - 1) * NumReps,
                drop = FALSE
            ], na.rm = TRUE)
        }
        if (profiles_scale) {
            MeanSet <- MeanSet - rowMeans(MeanSet, na.rm = TRUE)
        }

        layout(t(c(1, 1, 2, 2, 3, 3)))
        plot(0, 0,
            type = "n", bty = "n", xaxt = "n", yaxt = "n",
            xlab = NA, ylab = NA
        )
        legend("topright",
            col = rainbow(nrow(SubSet), alpha = 0.8, s = 0.7),
            legend = strtrim(rownames(SubSet), 20), lwd = 3,
            title = "Features"
        )
        plotCI(seq_len(NumCond) + runif(1, -0.1, 0.1), MeanSet[1, ],
            pch = 16,
            xlab = "Conditions", xlim = c(0.7, (NumCond + 1) - 0.7),
            ylab = "expression values",
            col = rainbow(nrow(MeanSet), alpha = 0.8, s = 0.7)[1],
            uiw = SDSet[1, ], type = "b", barcol = "#000000AA",
            ylim = range(cbind(MeanSet + SDSet, MeanSet - SDSet), na.rm = TRUE),
            xaxt = "none", lwd = 1.5
        )
        title(main = "Feature expression over conditions")
        axis(1, at = seq_len(NumCond), labels = colnames(MeanSet))
        # abline(h=fclim)
        if (nrow(MeanSet) > 1) {
            for (i in 2:nrow(MeanSet)) {
                plotCI(seq_len(NumCond) + runif(1, -0.1, 0.1),
                    MeanSet[i, ],
                    add = TRUE, pch = 16, col = rainbow(nrow(MeanSet))[i],
                    uiw = SDSet[i, ], type = "b", barcol = "#000000AA",
                    lwd = 1.5
                )
            }
        }

        if (length(SubSet) > 0) {
            message("Making circos plot")
            par(mar = rep(0, 4))
            circos.clear()
            circos.par(
                cell.padding = c(0, 0, 0, 0), canvas.xlim = c(-1.5, 1.5),
                canvas.ylim = c(-1.5, 1.5),
                track.margin = c(0, 0.02), start.degree = 90, gap.degree = 4
            )
            circos.initialize(seq_len(NumComps), xlim = c(0, 1))
            for (t in seq_len(NumTests)) {
                # print(tsign)
                nfeat <- min(nrow(SubSet), 30)
                cols <- rainbow(nfeat, alpha = 0.8, s = 0.7)
                tsign <- FCRegs[indices, (t - 1) * (NumComps) +
                    (seq_len(NumComps)),
                drop = FALSE
                ] < qlim
                circos.trackPlotRegion(
                    ylim = c(-3, 2), track.height = 1 / 12,
                    bg.border = "#777777",
                    panel.fun = function(x, y) {
                        name <- get.cell.meta.data("sector.index")
                        i <- get.cell.meta.data("sector.numeric.index")
                        xlim <- get.cell.meta.data("xlim")
                        ylim <- get.cell.meta.data("ylim")
                        xdiff <- (xlim[2] - xlim[1]) / nfeat
                        if (t == 1) {
                            circos.text(mean(xlim), max(ylim) + 30,
                                compNames[i],
                                facing = "inside",
                                niceFacing = TRUE, cex = 1, font = 2
                            )
                            circos.axis("top",
                                labels = strtrim(rownames(SubSetLR), 20),
                                major.at = seq(1 / (nfeat * 2), 1 - 1 /
                                    (nfeat * 2),
                                length = nfeat
                                ), minor.ticks = 0,
                                labels.cex = 0.8,
                                labels.facing = "reverse.clockwise",
                            )
                        }
                        for (j in which(tsign[, i])) {
                            circos.rect(
                                xleft = xlim[1] + (j - 1) * xdiff,
                                ybottom = ylim[1],
                                xright = xlim[2] - (nfeat - j) * xdiff,
                                ytop = ylim[2],
                                col = cols[j], border = NA
                            )
                        }
                    }
                )
            }
            fccols <- redblue(1001)
            # print(SubSet)
            circos.trackPlotRegion(
                ylim = c(-3, 2), track.height = 1 / 4,
                bg.border = NA, panel.fun = function(x, y) {
                    name <- get.cell.meta.data("sector.index")
                    i <- get.cell.meta.data("sector.numeric.index")
                    xlim <- get.cell.meta.data("xlim")
                    ylim <- get.cell.meta.data("ylim")
                    # circos.text(x=mean(xlim), y=1.7,
                    #            labels=name, facing = dd, cex=0.6,  adj = aa),
                    xdiff <- (xlim[2] - xlim[1]) / nfeat
                    for (j in seq_len(nfeat)) {
                        # print((SubSetLR[j,i]/max(LogRatios,na.rm=T))*500+500)
                        circos.rect(
                            xleft = xlim[1] + (j - 1) * xdiff,
                            ybottom = ylim[1],
                            xright = xlim[2] - (nfeat - j) * xdiff,
                            ytop = ylim[2],
                            col = fccols[(SubSetLR[j, i] /
                                max(LogRatios, na.rm = TRUE)) *
                                500 + 500], border = 0
                        )
                    }
                }
            )
            text(0, 0, "Log\nratios", cex = 0.7)
            # label the different tracks
            mtext(
                paste("Successful statistical tests\nfor threshold given above.
              \nFrom outer to inner circles",
                    paste("Track ", seq_len(NumTests), ": ", testNames,
                        sep = "",
                        collapse = "\n"
                    ),
                    sep = "\n"
                ),
                side = 1, outer = TRUE, adj = 1, line = -1, cex = 0.6
            )
            par(mar = c(5.1, 4.1, 4.1, .21))
        }
        message("plotting expression profiles finished")
    }
}

#' Plot UpSet
#'
#' Visualizes the intersections of significant features across multiple
#' comparisons using an UpSet plot. Summarizes all comparisons from all tests
#'
#' @param fulldata A `SummarizedExperiment` object containing the data.
#' @param qlim A numeric value, the q-value threshold for significance.
#' @param fclim A numeric vector, specifying fold change limits for filtering.
#'
#' @return An UpSet plot visualizing the intersections of significant features.
#' 
#' @examples
#' data(liver_example)
#'
#' plotUpset(liver_example, qlim = 0.05)
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom UpSetR upset
plotUpset <- function(fulldata, qlim = 0.05, fclim = c(0, 0)) {
    check_for_polystest(fulldata)

    # check compnames and testnames with metadata
    if (is.null(metadata(fulldata)$compNames) |
        is.null(metadata(fulldata)$testNames)) {
        stop("No metadata for statistical testing  found in the
        SummarizedExperiment object. Please run the statistical tests")
    }

    NumComps <- length(metadata(fulldata)$compNames)
    NumTests <- length(metadata(fulldata)$testNames)
    rdat <- rowData(fulldata)
    FCRegs <- filterFC(rdat, NumTests, NumComps, fclim)

    if (!is.null(FCRegs)) {
        WhereRegs <- FCRegs[, rep(
            0:(NumTests - 2),
            NumComps
        ) * (NumComps) + rep(seq_len(NumComps),
            each = NumTests - 1
        ),
        drop = FALSE
        ] < qlim
        WhereRegs[WhereRegs] <- 1
        deleted_cols <- which(colSums(WhereRegs, na.rm = TRUE) == 0)
        # print(deleted_cols)

        tcolnames <- paste("A", rep(seq_len(NumComps), each = NumTests - 1))
        if (length(deleted_cols) > 0) {
            tcolnames <- tcolnames[-deleted_cols]
            WhereRegs <- WhereRegs[, -deleted_cols, drop = FALSE]
        }
        tcols <- rep(rainbow(NumComps), each = 1)
        names(tcols) <- rep(paste("A", seq_len(NumComps)), 1)

        if (length(WhereRegs) > 0) {
            message("Plotting upset plots")
            upset_plot <- UpSetR::upset(as.data.frame(WhereRegs),
                nsets = ncol(WhereRegs),
                mainbar.y.label = "Significant features",
                nintersects = NA, keep.order = FALSE,
                sets = colnames(WhereRegs),
                text.scale = 1.5, mb.ratio = c(0.55, 0.45),
                set.metadata = list(
                    data = data.frame(
                        set = colnames(WhereRegs), cols = tcolnames,
                        crab = seq_len(ncol(WhereRegs))
                    ),
                    plots = list(list(
                        type = "matrix_rows", column = "cols",
                        colors = tcols, alpha = 0.5
                    ))
                )
            )
            message("Plotting upset plots finished")
            return(upset_plot)
        } else {
            NULL
        }
    }
}


#' Plot Number of Regulated Features
#'
#' @description This function plots the mfrnumber of regulated features across
#' comparisons for different statistical tests. It shows how the number of
#' significant features varies with different FDR thresholds.
#'
#' @param fulldata A `SummarizedExperiment` object containing the dataset.
#' @param compNames A character vector of comparison names. "all" selects all
#' comparisons.
#' @param testNames A character vector of test names used in the analysis.
#' Default values are "PolySTest", "limma", "Miss test", "rank products",
#' "permutation test", and "t-test".
#' @param qlim Numeric, q-value (FDR) threshold.
#' @param fclim Numeric vector, fold-change limits.
#' @param TestCols Character vector, colors to use for each test in the plot.
#' @param ... Arguments passed further to plot/lines calls
#'
#' @return Invisible. The function generates plots.
#'
#' @examples
#' data(liver_example)
#' plotRegNumber(fulldata = liver_example, NumComps = 3)
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData
#'
plotRegNumber <- function(fulldata, compNames = "all",
                          testNames = c(
                              "PolySTest", "limma", "Miss Test",
                              "rank products", "permutation test", "t-test"
                          ),
                          qlim = 0.05, fclim = c(0, 0),
                          TestCols = c(
                              "#33AAAA", "#33AA33", "#AA3333", "#AA33AA",
                              "#AAAA33", "#3333AA"
                          ), ...) {
    check_for_polystest(fulldata)
    rdat <- SummarizedExperiment::rowData(fulldata)

    # check compnames and testnames with metadata
    testNames2 <- make.names(testNames)
    testNames2 <- gsub("\\.", "_", testNames2)
    compNames <- check_stat_names(fulldata, compNames, testNames2)

    QValue <- as.matrix(rdat[, grep("^FDR_", colnames(rdat)), drop = FALSE])
    LogRatios <- as.matrix(rdat[, grep("^log_ratios_", colnames(rdat)),
        drop = FALSE
    ])
    # get only the p-values for the tests and comparisons
    QValue <- QValue[, grep(
        paste(paste0("^FDR_", testNames2), collapse = "|"),
        colnames(QValue)
    ), drop = FALSE]
    LogRatios <- LogRatios[, grep(paste(paste0("^log_ratios_", compNames),
        collapse = "|"
    ), colnames(LogRatios)),
    drop = FALSE
    ]
    QValue <- QValue[, grep(
        paste(paste0(compNames, "$"), collapse = "|"),
        colnames(QValue)
    ), drop = FALSE]

    NumTests <- length(testNames)
    NumComps <- length(compNames)

    FCRegs <- filterFC(cbind(LogRatios, QValue), NumTests, NumComps, fclim)

    if (!is.null(FCRegs)) {
        message("Plotting number of regulated features")
        set_mfrow(NumComps, 5)

        Qvalue <- as.matrix(rdat[, grep("^FDR", colnames(rdat)), drop = FALSE])
        tmpX <- 10^seq(log10(min(Qvalue, na.rm = TRUE)), 0.1, 0.01)
        tmpX[tmpX == 0] <- NA

        for (i in seq_len(NumComps)) {
            plotData <- vapply(
                tmpX, function(x) sum(FCRegs[, i] < x, na.rm = TRUE),
                numeric(1)
            )
            plotData[plotData == 0] <- NA

            plot(tmpX, plotData,
                main = paste("Comparison", i), xlab = "FDR threshold",
                ylab = "Number significant", type = "l",
                col = TestCols[1], ylim = c(1, nrow(Qvalue)), log = "xy",
                lwd = 2, ...
            )

            if (i == 1) {
                legend("topleft", legend = testNames, col = TestCols, lwd = 2)
            }

            for (j in seq_len(NumTests - 1)) {
                count_num <- vapply(
                    tmpX,
                    function(x) {
                        sum(FCRegs[, (NumComps) * j + i] < x,
                            na.rm = TRUE
                        )
                    }, numeric(1)
                )
                count_num[count_num == 0] <- NA
                lines(tmpX, count_num, col = TestCols[j + 1], lwd = 2, ...)
            }

            abline(v = qlim, col = "red")
        }

        message("Plotting number of regulated features finished")
    }
}

#' Heatmap Visualization with Heatmaply
#'
#' @description This function generates a heatmap for selected features across
#' comparisons using the heatmaply package.
#' It provides options for scaling and saving the plot to a file.
#'
#' @param fulldata A `SummarizedExperiment` object containing the dataset.
#' @param sel_prots Character vector specifying selected features to include in
#' the heatmap or "all" to include all proteins.
#' @param heatmap_scale Character, indicating if and how the data should be
#' scaled. Possible values are "none", "row", or "column".
#' @param file Optional character string specifying the path to save the heatmap
#' plot. If NULL, the plot is rendered interactively.
#' @param ... Arguments passed further to heatmaply function
#'
#' @return A plotly object if `file` is NULL. Otherwise, the heatmap is saved
#' to the specified file.
#'
#' @examples
#' data(liver_example)
#' plotHeatmaply(
#'     fulldata = liver_example, sel_prots = "all",
#'     heatmap_scale = "row"
#' )
#'
#' @importFrom plotly plotly_empty
#' @importFrom SummarizedExperiment rowData
#' @importFrom heatmaply heatmaply
#'
#' @export
plotHeatmaply <- function(fulldata, sel_prots = "all", heatmap_scale = "none",
                          file = NULL, ...) {
    check_for_polystest(fulldata)

    NumComps <- length(metadata(fulldata)$compNames)

    p <- plotly::plotly_empty()

    rdat <- SummarizedExperiment::rowData(fulldata)
    dat <- SummarizedExperiment::assay(fulldata)
    LogRatios <- as.matrix(rdat[, grep("^log_ratios_", colnames(rdat)),
        drop = FALSE
    ])
    Qvalue <- as.matrix(rdat[, grep("^FDR_", colnames(rdat)), drop = FALSE])
    rownames(LogRatios) <- rownames(Qvalue) <- rownames(rdat)

    if (ncol(Qvalue) > 0 & length(sel_prots) > 0) {
        if (all(sel_prots == "all")) {
            sel_prots <- seq_len(nrow(Qvalue))
        }

        NumCond <- metadata(fulldata)$NumCond
        NumReps <- metadata(fulldata)$NumReps

        # CI plots of max 30 features
        SubSetQval <- Qvalue[sel_prots, , drop = FALSE]
        SubSetLR <- LogRatios[sel_prots, , drop = FALSE]
        SubSetLR <- SubSetLR[
            order(rowMins(
                SubSetQval[, seq_len(NumComps),
                    drop = FALSE
                ],
                na.rm = TRUE
            )), ,
            drop = FALSE
        ]

        if (!is.null(SubSetLR)) {
            if (length(SubSetLR) > 0 & nrow(SubSetLR) > 1) {
                message("plotting heatmap ...")
                tdat <- dat[rownames(SubSetLR), (rep(
                    seq_len(NumReps),
                    NumCond
                ) - 1) *
                    NumCond + rep(seq_len(NumCond), each = NumReps),
                drop = FALSE
                ]
                rownames(tdat) <- strtrim(rownames(tdat), 30)

                # remove data rows with more than 45% missing values
                to_remove <- which(rowSums(is.na(tdat)) > ncol(tdat) * 0.45)
                tqvals <- Qvalue[rownames(SubSetLR), seq_len(NumComps),
                    drop = FALSE
                ]
                if (length(to_remove) > 0) {
                    tqvals <- tqvals[-to_remove, , drop = FALSE]
                    tdat <- tdat[-to_remove, , drop = FALSE]
                }
                # setting colors of p-values
                pcols <- rev(c(0.001, 0.01, 0.05, 1))
                ttt <- tqvals
                for (c in pcols) {
                    ttt[tqvals <= c] <- c
                }
                tqvals <- data.frame(ttt)
                for (c in seq_len(ncol(tqvals))) {
                    tqvals[, c] <- paste("<", as.character(tqvals[, c], pcols),
                        sep = ""
                    )
                }
                scaling <- heatmap_scale
                tqvals <- tqvals[order(rownames(tdat)), , drop = FALSE]
                p <- heatmaply::heatmaply(
                    tdat[order(rownames(tdat)), ,
                        drop = FALSE
                    ],
                    Colv = FALSE, scale = scaling, trace = "none", cexRow = 0.7,
                    plot_method = "plotly",
                    RowSideColors = tqvals, row_side_palette = grey.colors,
                    file = file, ...
                )
                message("Plotting heatmap finished")
            }
        }
    }
    p
}

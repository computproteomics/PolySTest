library(testthat)
library(SummarizedExperiment)

test_that("PolySTest_unpaired", {
  # Setup mock SummarizedExperiment object
  dataMatrix <- matrix(rnorm(1000), nrow = 100)
  colData <- DataFrame(Condition = rep(c("A", "B"), each = 5))
  rowData <- DataFrame(Gene = paste("Gene", 1:100))
  fulldata <- SummarizedExperiment(assay = list(quant = dataMatrix),
                                   colData = colData, rowData = rowData)
  metadata(fulldata) <- list(NumCond = 2, NumReps = 5)

  # Specifying comparisons
  allComps <- matrix(c("A", "B"), ncol = 2, byrow = TRUE)

  # Run function
  results <- PolySTest_unpaired(fulldata, allComps)

  # Check for added columns in rowData
  expect_true(any(grepl("p_values", colnames(rowData(results)))))
  expect_true(any(grepl("FDR", colnames(rowData(results)))))
  expect_equal(sum(rowData(results)$'FDR_PolySTest_B_vs_A' < 0.01), 0 )

  # Run with different setupt
  results <- PolySTest_unpaired(fulldata, allComps, statTests = c("t_test", "limma"))
  expect_true(any(grepl("FDR", colnames(rowData(results)))))
  expect_equal(sum(rowData(results)$'FDR_limma_B_vs_A' < 0.01, na.rm=T), 0 )

  # Run with different setupt
  results <- PolySTest_unpaired(fulldata, allComps, statTests = c("t_test", "limma", "rank_products"))
  expect_true(any(grepl("FDR", colnames(rowData(results)))))
  expect_equal(sum(rowData(results)$'FDR_PolySTest_B_vs_A' < 0.01), 0 )

  })

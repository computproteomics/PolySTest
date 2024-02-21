library(testthat)
library(PolySTest) # Replace with the name of your package

test_that("FindFCandQlim returns correct format", {
  # Mock data
  Qvalue <- matrix(runif(1800, 0, 0.05), ncol = 18)
  LogRatios <- matrix(rnorm(300, mean = 0, sd = 1), ncol = 3)

  # Expected result should be a numeric vector of length 2
  result <- FindFCandQlim(Qvalue, LogRatios)

  expect_type(result, "double")
  expect_length(result, 2)
})

test_that("FindFCandQlim returns correct format", {
  # Mock data
  data("liver_example")
  Qvalue <- as.data.frame(rowData(liver_example))
  LogRatios <- Qvalue[, grep("^log_ratios_",names(Qvalue))]
  Qvalue <- Qvalue[, grep("(^FDR_PolySTest_)|(^FDR_limma_)|(^FDR_Miss_Test_)",names(Qvalue))]

  # Expected result should be a numeric vector of length 2
  result <- FindFCandQlim(Qvalue, LogRatios)

  expect_type(result, "double")
  expect_length(result, 2)
})


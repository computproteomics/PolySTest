library(testthat)
library(PolySTest)

test_that("permtest_paired returns correct structure and types", {
  tMAData <- matrix(rnorm(100), nrow = 10)
  result <- permtest_paired(tMAData)

  expect_type(result, "list")
  expect_true(all(c("pPermutvalues", "qPermutvalues") %in% names(result)))
  expect_type(result$pPermutvalues, "double")
  expect_type(result$qPermutvalues, "double")
})


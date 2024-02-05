library(testthat)
# Assuming your package is named 'yourPackageName'
library(PolySTest)

context("Tests for Paired Function")

test_that("Paired function returns correct structure", {
  # Generating mock data
  MAData <- matrix(rnorm(600), ncol = 6) # 10 genes, 4 conditions
  NumCond <- 3
  NumReps <- 2

  # Running the function
  results <- Paired(MAData, NumCond, NumReps)

  # Testing return type
  expect_type(results, "list")

  # Testing return length
  expect_equal(length(results), 10)

  # Testing return structure
  expected_names <- c("lratios", "ptvalues", "plvalues", "pRPvalues", "pPermutvalues",
                      "qtvalues", "qlvalues", "qRPvalues", "qPermutvalues", "Sds")
  expect_equal(names(results), expected_names)

  # Testing content type of each component
  expect_type(results$lratios, "matrix")
  expect_type(results$ptvalues, "matrix")
  expect_type(results$plvalues, "matrix")
  # Add more as necessary for each component
})

# Add more tests here to cover edge cases, correctness of the analysis, etc.

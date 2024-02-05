library(testthat)
library(PolySTest) # Load your package

test_that("UnpairedDesign returns correct output structure", {
  # Creating mock data for the test
  Data <- matrix(rnorm(1000), ncol = 10) # 100 genes, 10 samples
  RRCateg <- matrix(c(1, 2, 2, 3), nrow = 2, ncol = 2) # Comparing conditions
  NumCond <- 5
  NumReps <- 2

  # Running the function
  result <- UnpairedDesign(Data, RRCateg, NumCond, NumReps)

  # Checking for correct output type
  expect_type(result, "list")

  # Checking for correct structure of the output
  expected_components <- c("lratios", "ptvalues", "plvalues", "pRPvalues",
                           "pPermutvalues", "qtvalues", "qlvalues",
                           "qRPvalues", "qPermutvalues", "Sds")
  expect_equal(names(result), expected_components)

  # Optionally, check the dimensions of the output components
  expect_equal(dim(result$lratios), c(nrow(Data), length(unique(RR[2, ]))-1))
  # Add more dimension checks for other components as appropriate
})

library(testthat)
context("Tests for Unpaired Function")

# Mock data and expected results for testing
test_data <- matrix(rnorm(600), ncol = 6)
num_cond <- 2
num_reps <- 3

# Test 1: Unpaired returns a list
test_that("Unpaired returns a list", {
  result <- Unpaired(test_data, num_cond, num_reps)
  expect_type(result, "list")
})

# Test 2: Unpaired returns list with correct components
test_that("Unpaired returns list with correct components", {
  result <- Unpaired(test_data, num_cond, num_reps)
  expected_names <- c("lratios", "ptvalues", "plvalues", "pRPvalues",
                      "pPermutvalues", "qtvalues", "qlvalues",
                      "qRPvalues", "qPermutvalues", "Sds")
  expect_equal(names(result), expected_names)
})

# Add more tests as necessary for each component

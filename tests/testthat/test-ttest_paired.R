library(testthat)
library(PolySTest)

test_that("ttest_paired returns correct structure and types", {
  tMAData <- matrix(rnorm(100), ncol = 5)
  result <- ttest_paired(tMAData)

  expect_true(all(c("ptvalues", "qtvalues") %in% names(result)))
  expect_type(result$ptvalues, "double")
  expect_type(result$qtvalues, "double")
})

test_that("ttest_paired handles NA values correctly", {
  tMAData_with_NA <- matrix(c(rnorm(1999), NA), ncol = 10)
  result <- ttest_paired(tMAData_with_NA)

  expect_true(all(is.na(result$ptvalues[nrow(tMAData_with_NA)])))
  expect_true(all(is.na(result$qtvalues[nrow(tMAData_with_NA)])))
})



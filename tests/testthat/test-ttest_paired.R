library(testthat)
library(PolySTest)

test_that("ttest_paired returns correct structure and types", {
  tMAData <- matrix(rnorm(1000), ncol = 5)
  result <- ttest_paired(tMAData)

  expect_true(all(c("ptvalues", "qtvalues") %in% names(result)))
  expect_type(result$ptvalues, "double")
  expect_type(result$qtvalues, "double")
})


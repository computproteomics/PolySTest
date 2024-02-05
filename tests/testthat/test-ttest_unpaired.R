test_that("ttest_unpaired returns correct output structure", {
  tData <- matrix(rnorm(100), nrow = 10)
  trefData <- matrix(rnorm(100), nrow = 10)
  result <- ttest_unpaired(tData, trefData)

  expect_type(result, "list")
  expect_true(all(c("ptvalues", "qtvalues") %in% names(result)))
  expect_type(result$ptvalues, "double")
  expect_type(result$qtvalues, "double")
})

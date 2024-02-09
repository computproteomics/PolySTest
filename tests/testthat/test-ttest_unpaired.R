test_that("ttest_unpaired returns correct output structure", {
  tData <- matrix(rnorm(1000), ncol = 10)
  trefData <- matrix(rnorm(1000), ncol = 10)
#  rownames(tData) <- rownames(trefData) <- paste("feature", 1:10)
  result <- ttest_unpaired(tData, trefData)

  expect_type(result, "list")
  expect_true(all(c("ptvalues", "qtvalues") %in% names(result)))
  expect_type(result$ptvalues, "double")
  expect_type(result$qtvalues, "double")
})

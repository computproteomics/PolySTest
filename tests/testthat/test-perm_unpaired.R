test_that("perm_unpaired returns correct output structure", {
  tData <- matrix(rnorm(10000), nrow = 100)
  trefData <- matrix(rnorm(10000), nrow = 100)
  result <- perm_unpaired(tData, trefData)

  expect_type(result, "list")
  expect_true(all(c("pPermutvalues", "qPermutvalues") %in% names(result)))
  expect_type(result$pPermutvalues, "double")
  expect_type(result$qPermutvalues, "double")
})

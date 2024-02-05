test_that("rp_unpaired returns correct output structure", {
  tData <- matrix(rnorm(1000), nrow = 100)
  trefData <- matrix(rnorm(1000), nrow = 100)
  result <- rp_unpaired(tData, trefData)

  expect_type(result, "list")
  expect_true(all(c("pRPvalues", "qRPvalues") %in% names(result)))
  expect_type(result$pRPvalues, "double")
  expect_type(result$qRPvalues, "double")
})

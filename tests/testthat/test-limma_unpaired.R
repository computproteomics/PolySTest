library(testthat)
library(PolySTest)

test_that("limma_unpaired returns correct output structure", {
  data <- matrix(rnorm(1500), ncol = 15)
  RRCateg <- matrix(c(1, 2, 2, 3), nrow = 2, ncol = 2)
  result <- limma_unpaired(data, 3, 5, RRCateg)

  expect_type(result, "list")
  expect_true(all(c("plvalues", "qlvalues", "Sds") %in% names(result)))
  expect_true(is.matrix(result$plvalues))
  expect_true(is.matrix(result$qlvalues))
  expect_true(is.vector(result$Sds))
})

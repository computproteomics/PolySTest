library(SummarizedExperiment)
test_that("plotRegNumber runs without error", {
  data(liver_example)
  # Remove margins to avoid figure margins too large error
  expect_message(plotRegNumber(fulldata = liver_example,
                                 NumComps = 3),
                 regexp="finished")
})

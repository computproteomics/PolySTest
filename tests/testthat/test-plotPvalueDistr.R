test_that("plotPvalueDistr runs without error", {
  data(liver_example)
  # Remove margins to avoid figure margins too large error
  expect_message(plotPvalueDistr(fulldata = liver_example,
                                 compNames = c("FO.Rep._vs_HF.Rep", "FO.Rep._vs_HF.Rep",
                                               "TTA.FO.Rep._vs_HF.Rep."),
                                 testCols = c("#33AAAA")),
                 regexp="finished")
  })

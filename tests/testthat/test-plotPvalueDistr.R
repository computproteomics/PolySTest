test_that("plotPvalueDistr runs without error", {
  data(liver_example)
  # Remove margins to avoid figure margins too large error
  expect_message(plotPvalueDistr(fulldata = liver_example, c("Miss Test","t-test"),
                                 compNames = c("HF.Rep._vs_TTA.Rep.", "FO.Rep._vs_TTA.Rep.",
                                               "TTA.FO.Rep._vs_TTA.Rep."),
                                               testCols = c("#33AAAA")),
                 regexp="finished")
  })

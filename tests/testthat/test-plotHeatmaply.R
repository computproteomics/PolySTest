test_that("plotHeatmaply handles inputs correctly", {
  data("liver_example")
  # Test: Call plotVolcano with the mock data and check for errors or warnings
  expect_message(plotHeatmaply(liver_example,
                             3,
                             sel_prots = "all"),
                 regexp = "finished")
  expect_message(plotHeatmaply(liver_example,
                             3,
                             sel_prots = 5:20,
                             heatmap_scale = "row"),
                 regexp = "finished")
})

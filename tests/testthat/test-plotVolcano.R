test_that("plotVolcano handles inputs correctly", {
  data("liver_example")
  # Test: Call plotVolcano with the mock data and check for errors or warnings
  expect_message(plotVolcano(liver_example,
                             compNames = c("HF.Rep._vs_TTA.Rep.", "FO.Rep._vs_TTA.Rep.",
                                           "TTA.FO.Rep._vs_TTA.Rep."
                             ),
                             sel_prots = "all",
                             qlim = 0.05, fclim = c(-2, 2)),
                 regexp = "finished")
  expect_message(plotVolcano(liver_example,
                             compNames = c("HF.Rep._vs_TTA.Rep.", "FO.Rep._vs_TTA.Rep.",
                                           "TTA.FO.Rep._vs_TTA.Rep."
                             ),
                             sel_prots = 5:20,
                             qlim = 0.05, fclim = c(-2, 2)),
                 regexp = "finished")
})

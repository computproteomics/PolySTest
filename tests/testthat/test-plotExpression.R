test_that("plotExpression executes without error", {
  data("liver_example")
  # Test: Call plotVolcano with the mock data and check for errors or warnings
  expect_message(suppressWarnings(plotExpression(liver_example,
                                                 compNames = c("HF.Rep._vs_TTA.Rep.", "FO.Rep._vs_TTA.Rep.",
                                                               "TTA.FO.Rep._vs_TTA.Rep."
                                                 ),
                             sel_prots = 1:10,
                             qlim = 0.05, fclim = c(-1, 1))),
                 regexp = "finished")
  expect_message(suppressWarnings(plotExpression(liver_example,
                                                 compNames = c("HF.Rep._vs_TTA.Rep.", "FO.Rep._vs_TTA.Rep.",
                                                               "TTA.FO.Rep._vs_TTA.Rep."
                                                 ),
                             sel_prots = 5:20,
                            profiles_scale = FALSE,
                             qlim = 0.05, fclim = c(-2, 2))),
                 regexp = "finished")

  })

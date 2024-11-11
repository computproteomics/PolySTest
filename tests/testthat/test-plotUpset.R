test_that("plotUpset executes without error", {
    data("liver_example")
  # Test: Call plotVolcano with the mock data and check for errors or warnings
  expect_message(plotUpset(liver_example, qlim = 0.05, fclim = c(-1, 1)),
                 regexp = "finished")
})

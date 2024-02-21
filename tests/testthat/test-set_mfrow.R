test_that("set_mfrow sets correct layout", {
  expect_silent(set_mfrow(num_total = 6, max_col = 3))
  layout_matrix <- par("mfrow")
  expect_equal(layout_matrix, c(2, 3))
})

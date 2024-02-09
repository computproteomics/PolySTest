library(testthat)

test_that("update_conditions_with_lcs updates condition names correctly", {
  # Assuming 'se' and 'default_conditions' are defined as above
   se <- SummarizedExperiment(assays = list(count = matrix(rnorm(200), ncol = 10)))
  metadata(se) <- list(NumCond = 2, NumReps = 5)
  rownames(colData(se)) <-  paste0(rep(c("CondA_Rep","CondB_Rep"), 5),
                                    rep(1:5, each=2))
  default_conditions <- c("Condition_A", "Condition_B")
  expected_conditions <- c("CondA_Rep", "CondB_Rep") # Expected output based on mock data
  updated_conditions <- update_conditions_with_lcs(se, default_conditions)

  expect_equal(unique(colData(updated_conditions)$Condition), expected_conditions)
})

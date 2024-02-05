test_that("overall test of full workflow", {
  # package raw data file
  filename <- system.file("extdata", "LiverAllProteins.csv", package = "PolySTest")
  NumReps <- 3
  NumCond <- 4
  isPaired <- TRUE

  ## Read file
  dat <- read.csv(filename, row.names=1, stringsAsFactors = F)
  ndatcol <- ncol(dat)


  # Normalize
  dat <- t(t(dat) -
             colMedians(as.matrix(dat), na.rm = TRUE))

  ####### Defining conditions from names
  conditions <- update_conditions_with_lcs(dat, NumCond, NumReps, conditions = paste("C", 1:NumCond, sep = ""))

  ####### Defining comparison for statistical tests
  FullReg <- allComps <- NULL

  allComps <- create_pairwise_comparisons(conditions, 1)
  ncomps <- nrow(allComps)

  RR <- convert_comps_to_indices(allComps, conditions, NumCond, NumReps)

  # for maybe later to include onlyLIMMA option
  MAData <- NULL
  for (i in seq_len(ncol(RR))) {
    MAData <- cbind(MAData, dat[, RR[1, i]] - dat[, RR[2, i]])
  }
  rownames(MAData) <- rownames(dat)
  qvalues <- Paired(MAData, ncomps, NumReps)
  qvalues <- UnpairedDesign(dat, RR, NumCond, NumReps)
  MissingStats <- list(
    pNAvalues = matrix(1, ncol = ncomps, nrow = nrow(dat)),
    qNAvalues = matrix(1, ncol = ncomps, nrow = nrow(dat))
  )
  MissingStats <- MissingStatsDesign(dat, RR, NumCond, NumReps)

  # Prepare output data
  FullReg <- prepare_output_data(dat, qvalues, MissingStats, allComps)
  DRFs <- FullReg[,32] < 0.01
  expect_equal(sum(DRFs), 21)
  expect_match(names(DRFs)[4], "A7VJC2")

})


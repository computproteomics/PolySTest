test_that("overall test of full workflow", {
  library(SummarizedExperiment)
  # package raw data file
  filename <- system.file("extdata", "LiverAllProteins.csv", package = "PolySTest")
  NumReps <- 3
  NumCond <- 4
  isPaired <- TRUE

  ## Read file
  dat <- read.csv(filename, row.names=1, stringsAsFactors = F)

  ## Reduce to first 200
  dat <- dat[1:200,]

  sampleMetadata <- data.frame(Condition = rep(paste("Condition", 1:NumCond), NumReps),
                               Replicate = rep(1:NumReps, each=NumCond))

  # Create the SummarizedExperiment object
  fulldata <- SummarizedExperiment(assays = list(quant = dat),
                                   colData = sampleMetadata)
  # Adding metadata
  metadata(fulldata) <- list(
    NumReps = NumReps,
    NumCond = NumCond
  )

  # Access the assay data
  dat <- assay(fulldata, "quant")

  # Normalize
  dat <- t(t(dat) -
             colMedians(as.matrix(dat), na.rm = TRUE))

  # Update the assay data in 'fulldata'
  assay(fulldata, "quant") <- dat

  ####### Defining conditions from names
  fulldata <- PolySTest:::update_conditions_with_lcs(fulldata)
  conditions <- unique(colData(fulldata)$Condition)

  ####### Defining comparison for statistical tests
  FullReg <- allComps <- NULL

  allComps <- PolySTest:::create_pairwise_comparisons(conditions, 2)

  # Run paired tests
  fulldata <- PolySTest_paired(fulldata, allComps)

  DRFs <- rowData(fulldata)[,"FDR_PolySTest_HF.Rep._vs_TTA.Rep."] < 0.01
  expect_equal(sum(DRFs), 23)
  expect_match(rownames(rowData(fulldata))[4], "A7VJC2")

  # Run unpaired tests
  fulldata <- PolySTest_unpaired(fulldata, allComps)
  DRFs <- rowData(fulldata)[,"FDR_PolySTest_HF.Rep._vs_TTA.Rep."] < 0.01
  expect_equal(sum(DRFs), 23)
  expect_match(rownames(rowData(fulldata))[4], "A7VJC2")


})


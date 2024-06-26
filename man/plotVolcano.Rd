% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualizations.R
\name{plotVolcano}
\alias{plotVolcano}
\title{Plot Volcano Plots for PolySTest results}
\usage{
plotVolcano(
  fulldata,
  compNames = "all",
  testNames = c("PolySTest", "limma", "Miss Test", "rank products", "permutation test",
    "t-test"),
  sel_prots = "all",
  qlim = 0.05,
  fclim = c(0, 0),
  testCols = c("#33AAAA", "#33AA33", "#AA3333", "#AA33AA", "#AAAA33", "#3333AA"),
  ...
)
}
\arguments{
\item{fulldata}{A \code{SummarizedExperiment} object containing the data.}

\item{compNames}{A character vector of comparison names.}

\item{testNames}{A character vector of test names including "PolySTest",
"limma", "Miss test", "rank products", "permutation test", and "t-test".}

\item{sel_prots}{A numeric vector indicating selected features to be
visualized differently or "all" to select all features. Default is "all".}

\item{qlim}{A numeric value setting the q-value limit for the plots.
Default is 0.05.}

\item{fclim}{A numeric vector of length two setting the fold-change limits
for the plots. Default is c(0,0).}

\item{testCols}{A character vector of colors for each test. Default is a
predefined set of colors.}

\item{...}{Additional arguments passed to the plot function.}
}
\value{
Creates volcano plots for the specified tests and comparisons.
}
\description{
This function creates volcano plots for all specified statistical tests and
comparisons using data from a \code{SummarizedExperiment} object. It highlights
selected proteins and applies fold-change and q-value limits for
visualization.
}
\examples{
data(liver_example)
compNames <- c("HF.Rep._vs_TTA.Rep.")
plotVolcano(liver_example, compNames)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualizations.R
\name{plotRegNumber}
\alias{plotRegNumber}
\title{Plot Number of Regulated Features}
\usage{
plotRegNumber(
  fulldata,
  compNames = "all",
  testNames = c("PolySTest", "limma", "Miss Test", "rank products", "permutation test",
    "t-test"),
  qlim = 0.05,
  fclim = c(0, 0),
  TestCols = c("#33AAAA", "#33AA33", "#AA3333", "#AA33AA", "#AAAA33", "#3333AA"),
  ...
)
}
\arguments{
\item{fulldata}{A \code{SummarizedExperiment} object containing the dataset.}

\item{compNames}{A character vector of comparison names. "all" selects all
comparisons.}

\item{testNames}{A character vector of test names used in the analysis.
Default values are "PolySTest", "limma", "Miss test", "rank products",
"permutation test", and "t-test".}

\item{qlim}{Numeric, q-value (FDR) threshold.}

\item{fclim}{Numeric vector, fold-change limits.}

\item{TestCols}{Character vector, colors to use for each test in the plot.}

\item{...}{Arguments passed further to plot/lines calls}
}
\value{
Invisible. The function generates plots.
}
\description{
This function plots the mfrnumber of regulated features across
comparisons for different statistical tests. It shows how the number of
significant features varies with different FDR thresholds.
}
\examples{
data(liver_example)
plotRegNumber(fulldata = liver_example, NumComps = 3)

}

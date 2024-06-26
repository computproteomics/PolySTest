% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFuncs.R
\name{FindFCandQlimAlternative}
\alias{FindFCandQlimAlternative}
\title{Function to determine "optimal" fold-change and q-value thresholds using
a higher criticism method and an FC threshold based on a standard deviation
of one}
\usage{
FindFCandQlimAlternative(Pvalue, LogRatios)
}
\arguments{
\item{Pvalue}{A matrix of p-values for each condition}

\item{LogRatios}{A matrix of log ratios for each condition}
}
\value{
A vector containing the mean fold-change threshold and mean q-value
threshold
}
\description{
Function to determine "optimal" fold-change and q-value thresholds using
a higher criticism method and an FC threshold based on a standard deviation
of one
}
\examples{
Pvalue <- matrix(seq(0.01, 0.12, 0.01), nrow = 4, ncol = 3)
LogRatios <- matrix(c(
    1.2, 0.8, 1.5, -0.5, 0.2, 0.9, -1.1, 0.7, 1.8,
    -0.9, 0.3, 1.1
), nrow = 4, ncol = 3)
thresholds <- FindFCandQlimAlternative(Pvalue, LogRatios)
print(thresholds)

}

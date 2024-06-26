% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFuncs.R
\name{FindFCandQlim}
\alias{FindFCandQlim}
\title{Find Optimal Fold-Change and Q-Value Thresholds}
\usage{
FindFCandQlim(Qvalue, LogRatios)
}
\arguments{
\item{Qvalue}{A matrix of q-values obtained from statistical tests such as
PolySTest, limma, rank products, and permutation tests. NA values in the
matrix are treated as non-significant. The matrix should have the same number
of rows as the LogRatios matrix, and the number of columns should be a
multiple of the number of conditions, given by the number of different
statistical tests.}

\item{LogRatios}{A matrix of log2 fold-change ratios for the same
comparisons.}
}
\value{
A numeric vector containing two elements: the optimal fold-change
threshold and the optimal q-value threshold, which together maximize the
number of significant features detected.
}
\description{
Identifies the optimal fold-change (FC) and q-value thresholds that maximize
the number of significant features identified across different statistical
tests. It processes a matrix of q-values, applying various fold-change
thresholds,
and computes the distribution of significant features for each q-value
threshold
to determine the optimal combination of thresholds.
}
\details{
The function first replaces NA values in the Qvalue matrix with 1 to
denote non-significant results. It then identifies the smallest q-value
and calculates a range of q-values around this minimum. For each fold-change
threshold, the function adjusts the Qvalue matrix and calculates the number
of significant features for each q-value threshold. The combination yielding
the highest mean distribution of significant features is considered optimal.
}
\examples{
# Example Qvalue and LogRatios matrices
Qvalue <- matrix(seq(0.01, 0.12, 0.01), nrow = 4, ncol = 3)
LogRatios <- matrix(c(
    1.2, 0.8, 1.5, -0.5, 0.2, 0.9, -1.1, 0.7,
    1.8, -0.9, 0.3, 1.1
), nrow = 4, ncol = 3)
# Find optimal thresholds
thresholds <- FindFCandQlim(Qvalue, LogRatios)
print(thresholds)

}
\keyword{internal}

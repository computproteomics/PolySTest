% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualizations.R
\name{set_mfrow}
\alias{set_mfrow}
\title{Set Graphics Layout}
\usage{
set_mfrow(num_total, max_col)
}
\arguments{
\item{num_total}{The total number of plots to display.}

\item{max_col}{The maximum number of columns for the layout.}
}
\value{
Sets the \code{mfrow} parameter for \code{par()}.
}
\description{
Automatically sets the \code{mfrow} parameter for \code{par()} based on the total
number of plots and the maximum number of columns desired.
}
\examples{
# This will set the layout to 2 rows of 3 columns
set_mfrow(num_total = 6, max_col = 3)
for (i in 1:6) hist(1:10)
par(mfrow = c(1, 1))
}

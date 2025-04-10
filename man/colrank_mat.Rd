% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_functions.R
\name{colrank_mat}
\alias{colrank_mat}
\title{Rank the columns of a matrix such that rank=1 is the highest positive value}
\usage{
colrank_mat(mat, ties_arg = "min", na_arg = "keep")
}
\arguments{
\item{mat}{A dense numeric matrix}

\item{ties_arg}{A character indicating how ties are assigned, as fed into
the ties argument of the base rank() function. "min" is the assumed default.}

\item{na_arg}{A character indicating how NAs are handled, as fed into
the na.last argument of the base rank() function. "keep" is the assumed default.}
}
\value{
Return a matrix of integer ranks the same dimension as input matrix
}
\description{
Rank the columns of a matrix such that rank=1 is the highest positive value
}

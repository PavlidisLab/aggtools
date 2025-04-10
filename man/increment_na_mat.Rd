% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_functions.R
\name{increment_na_mat}
\alias{increment_na_mat}
\title{Count the NAs in cmat and increment the corresponding indices of na_mat}
\usage{
increment_na_mat(cmat, na_mat)
}
\arguments{
\item{cmat}{A dense gene by gene correlation matrix}

\item{na_mat}{A dense gene by gene matrix that tracks NA presence across
multiple correlation matrices}
}
\value{
na_mat with the indices corresponding to NAs in cmat increased by 1
}
\description{
Count the NAs in cmat and increment the corresponding indices of na_mat
}

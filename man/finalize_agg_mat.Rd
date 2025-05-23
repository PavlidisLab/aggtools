% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_functions.R
\name{finalize_agg_mat}
\alias{finalize_agg_mat}
\title{Finalize aggregate matrix}
\usage{
finalize_agg_mat(amat, agg_method, n_celltypes, na_mat)
}
\arguments{
\item{amat}{A dense gene by gene matrix tracking the aggregate correlations}

\item{agg_method}{One of "allrank", "colrank", or "FZ"}

\item{n_celltypes}{The count of cell types that went into the aggregated
correlation}

\item{na_mat}{A dense gene by gene matrix tracking the count of NAs across
the correlation matrices}
}
\value{
A dense gene by gene matrix of the formatted aggregated correlations
}
\description{
Format the final summed aggregate correlations depending on the aggregate
strategy.
}
\details{
all_rank: jointly re-rank the summed ranks and standardize back into
[0, 1], then convert back to a symmetric matrix for downstream operations

col_rank: re-rank the summed ranks for each column separately and standardize
back into [0, 1]

FZ: divide each element by its count of its measured (i.e., non-NA)
observations
}

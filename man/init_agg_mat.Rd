% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_functions.R
\name{init_agg_mat}
\alias{init_agg_mat}
\title{Initiate aggregate matrix}
\usage{
init_agg_mat(pc_df)
}
\arguments{
\item{pc_df}{A data frame of unique protein coding gene symbols}
}
\value{
A gene by gene matrix of 0s with dimensions equal to the count of
symbols in pc_df
}
\description{
Initates a gene by gene matrix of 0s for holding aggregate correlations.
pc_df assumed to be protein coding table with "Symbol" as a column and unique
genes within Symbol.
}

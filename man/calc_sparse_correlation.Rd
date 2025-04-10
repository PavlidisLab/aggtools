% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_functions.R
\name{calc_sparse_correlation}
\alias{calc_sparse_correlation}
\title{Wrapper for calling sparse Pearson or Spearman correlation}
\usage{
calc_sparse_correlation(mat, cor_method)
}
\arguments{
\item{mat}{A sparse numeric m by n Matrix}

\item{cor_method}{One of "pearson" or "spearman"}
}
\value{
A dense n by n matrix of the correlations between the columns of mat
}
\description{
Wrapper for calling sparse Pearson or Spearman correlation
}

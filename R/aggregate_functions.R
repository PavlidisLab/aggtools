#' Rank the columns of a matrix such that rank=1 is the highest positive value
#'
#' @param mat A numeric matrix
#' @param ties_arg A character indicating how ties are assigned, as fed into
#' the ties argument of the base rank() function. "min" is the assumed default.
#' @param na_arg A character indicating how NAs are handled, as fed into
#' the na.last argument of the base rank() function. "keep" is the assumed default.
#'
#' @return Return a matrix of integer ranks the same dimension as input matrix
#' @export
#'
#' @examples
colrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {

  rank_mat <- apply(-mat, 2, rank, ties.method = ties_arg, na.last = na_arg)
  return(rank_mat)
}




#' Rank the entire matrix jointly such that rank=1 is the highest positive value
#'
#' @param mat A numeric matrix
#' @param ties_arg A character indicating how ties are assigned, as fed into
#' the ties argument of the base rank() function. "min" is the assumed default.
#' @param na_arg A character indicating how NAs are handled, as fed into
#' the na.last argument of the base rank() function. "keep" is the assumed default.
#'
#' @return Return a matrix of integer ranks the same dimension as input matrix
#' @export
#'
#' @examples
allrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {

  rank_mat <- rank(-mat, ties.method = ties_arg, na.last = na_arg)
  rank_mat <- matrix(rank_mat, nrow = nrow(mat), ncol = ncol(mat))
  rownames(rank_mat) <- colnames(rank_mat) <- rownames(mat)
  return(rank_mat)
}




# TODO: consider trying to force to sparse?
# Wraps to qlcMatrix::sparseCor() to make resulting matrix inherit input names

#' Column-wise Pearson's correlation for sparse matrices
#'
#' @importFrom qlcMatrix corSparse
#' @import Matrix
#' @param mat A sparse numeric m by n Matrix
#' @return A dense n by n matrix of the Pearson's correlation between the
#' columns of mat
#' @export
#'
#' @examples
sparse_pcor <- function(mat) {

  stopifnot(inherits(mat, "dgCMatrix"))

  cmat <- qlcMatrix::corSparse(mat)
  colnames(cmat) <- rownames(cmat) <- colnames(mat)

  return(cmat)
}

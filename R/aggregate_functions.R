#' Rank the columns of a matrix
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

# TODO: consider trying to force to sparse in cor call?


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



#' Wrapper for calling sparse Pearson or Spearman correlation
#'
#' @param mat A sparse numeric m by n Matrix
#' @param cor_method One of "pearson" or "spearman"
#'
#' @return A dense n by n matrix of the correlations between the columns of mat
#' @export
#'
#' @examples
calc_sparse_correlation <- function(mat, cor_method) {

  stopifnot(inherits(mat, "dgCMatrix"))
  stopifnot(cor_method %in% c("pearson", "spearman"))

  cmat <- if (cor_method == "pearson") {
    sparse_pcor(mat)
  } else {
    sparse_scor(mat)
  }

  return(cmat)
}



#' Set NAs in matrix to 0
#'
#' @param mat A matrix
#'
#' @return A matrix of the same dimensions where all values are preserved save
#' for NAs assigned to 0
#' @export
#'
#' @examples
na_to_zero <- function(mat) {

  stopifnot(is.matrix(mat))
  mat[is.na(mat)] <- 0
  return(mat)
}



#' Set diag of a square matrix to 1
#'
#' @param mat A matrix
#'
#' @return A matrix of the same dimensions where all values are preserved save
#' for the diagonals assigned to 1
#' @export
#'
#' @examples
diag_to_one <- function(mat) {

  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  diag(mat) <- 1
  return(mat)
}



#' Set upper triangle of a square matrix to NA
#'
#' @param mat A matrix
#'
#' @return A matrix of the same dimensions where all values are preserved save
#' for the values above the diagonal set to NA. Diagonal is preserved.
#' @export
#'
#' @examples
uppertri_to_na <- function(mat) {

  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  mat[upper.tri(mat)] <- NA
  return(mat)
}



#' Replace the upper triangle of square matrix with its lower triangle.
#'
#' @param mat A matrix
#'
#' @return A matrix of the same dimensions where all upper triangular values are
#' replaced by its lower triangle values.
#' @export
#'
#' @examples
lowertri_to_symm <- function(mat) {

  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  mat[upper.tri(mat)] <-  t(mat)[upper.tri(mat)]
  return(mat)
}



#' Set sparse columns to 0.
#'
#'# If a col of mat has fewer non-zero elements than min_count, set that col to
#'0. This is done to produce an NA during correlation, instead of allowing cors
#' derived from overly sparse columns.
#' n=20 default used from https://doi.org/10.1093/bioinformatics/btv118
#'
#' @import Matrix
#' @param mat A sparse matrix
#' @param min_count A non-negative integer. 20 is the assumed default.
#'
#' @return
#' @export
#'
#' @examples
zero_sparse_cols <- function(mat, min_count = 20) {

  stopifnot(inherits(mat, "dgCMatrix"),
            is.numeric(min_count),
            min_count >= 0 & min_count <= nrow(mat))

  nonzero_cells <- colSums(mat != 0)
  filt_genes <- nonzero_cells < min_count
  if (any(filt_genes)) mat[, filt_genes] <- 0

  return(mat)
}

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


#' Fisher Z transformation
#'
#' @param cor A numeric correlation value (single element, array, matrix)
#'
#' @return The transformed correlations
#' @seealso DescTools::FisherZ()
#' @export
#'
#' @examples
fisherz <- function(cor) {

  stopifnot(is.numeric(cor), all(cor >= -1 & cor <= 1))

  fz <- 0.5 * log((1 + cor) / (1 - cor))
  return(fz)
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
#' for the values above the diagonal set to NA and the diagonal is preserved
#' @export
#'
#' @examples
uppertri_to_na <- function(mat) {

  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  mat[upper.tri(mat)] <- NA
  return(mat)
}



#' Replace the upper triangle of square matrix with its lower triangle
#'
#' @param mat A matrix
#'
#' @return A matrix of the same dimensions where all upper triangular values are
#' replaced by its lower triangle values
#' @export
#'
#' @examples
lowertri_to_symm <- function(mat) {

  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  mat[upper.tri(mat)] <-  t(mat)[upper.tri(mat)]
  return(mat)
}



#' Set sparse columns to 0
#'
#'# If a col of mat has fewer non-zero elements than min_count, set that col to
#'0. This is done to produce an NA during correlation, instead of allowing cors
#' derived from overly sparse columns.
#' n=20 default used from https://doi.org/10.1093/bioinformatics/btv118
#'
#' @import Matrix
#' @param mat A sparse matrix
#' @param min_count A non-negative integer - 20 is the assumed default
#'
#' @return A matrix of the same dimensions where all sparse columns are set to 0
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



#' Prepare a cell type matrix for coexpression
#'
#' Subset mat to cells annotated to the input cell_type, set sparse genes to 0
#' and transpose the resulting matrix.

#' @param mat A sparse gene by cell matrix
#' @param meta A data frame that maps cell IDs to cell types
#' @param cell_type A character of the cell type to subset mat on
#' @param min_count A non-negative integer - 20 is the assumed default
#'
#' @return A sparse cell by gene matrix
#' @export
#'
#' @examples
prepare_celltype_mat <- function(mat, meta, cell_type, min_count = 20) {

  stopifnot(inherits(mat, "dgCMatrix"),
            c("ID", "Cell_type") %in% colnames(meta),
            cell_type %in% meta[["Cell_type"]])

  ids <- meta[meta$Cell_type %in% cell_type, "ID"]
  stopifnot(all(ids %in% colnames(mat)))

  ct_mat <- t(mat[, ids])
  ct_mat <- zero_sparse_cols(ct_mat, min_count)
  stopifnot(all(rownames(ct_mat) %in% meta[["ID"]]))

  return(ct_mat)
}



#' Initiate aggregate matrix
#'
#' Initates a gene by gene matrix of 0s for holding aggregate correlations.
#' pc_df assumed to be protein coding table with "Symbol" as a column and unique
#' genes within Symbol.
#'
#' @param pc_df A data frame of unique protein coding gene symbols
#'
#' @return A gene by gene matrix of 0s with dimensions equal to the count of
#' symbols in pc_df
#' @export
#'
#' @examples
init_agg_mat <- function(pc_df) {

  stopifnot("Symbol" %in% colnames(pc_df),
            length(unique(pc_df[["Symbol"]])) == length(pc_df[["Symbol"]]))

  amat <- matrix(0, nrow = nrow(pc_df), ncol = nrow(pc_df))
  rownames(amat) <- colnames(amat) <- pc_df[["Symbol"]]
  return(amat)
}



#' Count the NAs in cmat and increment the corresponding indices of na_mat
#'
#' @param cmat A gene by gene correlation matrix
#' @param na_mat A gene by gene matrix that tracks NA presence across multiple
#' correlation matrices
#'
#' @return na_mat with the indices corresponding to NAs in cmat increased by 1
#' @export
#'
#' @examples
increment_na_mat <- function(cmat, na_mat) {

  na_ix <- which(is.na(cmat), arr.ind = TRUE)
  na_mat[na_ix] <- na_mat[na_ix] + 1
  return(na_mat)
}



#' Transform correlation matrix
#'
#' Take an input correlation matrix and transform it to prepare for aggregation.
#' Transformation entails all ranking, column ranking, or Fisher's Z.
#' Each methods impute NA cors to 0 and ensure that the diagonal equals 1.
#'
#' allrank: make the matrix tri. to prevent double ranking symmetric elements
#' and then jointly rank the matrix (lower rank = positive cor).
#'
#' colrank: rank each column separately (lower rank = positive cor).
#'
#' FZ: perform Fisher's Z transform on the raw correlations

#' @param cmat A gene by gene correlation matrix
#' @param agg_method One of "allrank", "colrank", or "FZ"
#'
#' @return A gene by gene matrix of transformed correlations
#' @export
#'
#' @examples
transform_correlation_mat <- function(cmat, agg_method) {

  stopifnot(agg_method %in% c("allrank", "colrank", "FZ"))
  stopifnot(is.matrix(cmat), identical(rownames(cmat), colnames(cmat)))

  cmat <- cmat %>%
    na_to_zero() %>%
    diag_to_one()

  if (agg_method == "allrank") {

    cmat <- cmat %>%
      uppertri_to_na() %>%
      allrank_mat()

  } else if (agg_method == "colrank") {

    cmat <- colrank_mat(cmat)

  } else if (agg_method == "FZ") {

    cmat <- fisherz(cmat)
  }

  return(cmat)
}



#' Finalize aggregate matrix
#'
#' Format the final summed aggregate correlations depending on the aggregate
#' strategy.
#'
#' all_rank: jointly re-rank the summed ranks and standardize into [0, 1], then
#' convert back to a symmetric matrix for ease of downstream operations
#'
#' col_rank: re-rank the summed ranks for each column separately and standardize
#' into [0, 1]
#'
#' FZ: divide each element by its count of its measured (i.e., non-NA)
#' observations
#'
#' @param amat A gene by gene matrix tracking the aggregate correlations
#' @param agg_method One of "allrank", "colrank", or "FZ"
#' @param n_celltypes The count of cell types that went into the aggregated
#' correlation
#' @param na_mat A gene by gene matrix tracking the count of NAs across the
#' correlation matrices
#'
#' @return A gene by gene matrix of the formatted aggregated correlations
#' @export
#'
#' @examples
finalize_agg_mat <- function(amat, agg_method, n_celltypes, na_mat) {

  stopifnot(agg_method %in% c("allrank", "colrank", "FZ"))
  stopifnot(is.matrix(amat), identical(rownames(amat), colnames(amat)))

  if (agg_method == "allrank") {

    amat <- allrank_mat(amat) / sum(!is.na(amat))
    amat <- diag_to_one(amat)
    amat <- lowertri_to_symm(amat)

  } else if (agg_method == "colrank") {

    amat <- colrank_mat(amat)
    ngene <- nrow(amat)
    amat <- apply(amat, 2, function(x) x/ngene)

  } else if (agg_method == "FZ") {

    amat <- amat / (n_celltypes - na_mat)

  }

  return(amat)
}

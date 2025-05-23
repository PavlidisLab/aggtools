# TODO: more explicit check of NA content in sparse_cor()
# TODO: revisit aggr_coexpr_single_dataset test of insufficient counts


# Helper to create mock data
# ------------------------------------------------------------------------------



generate_test_data <- function() {

  set.seed(123)

  # 100 cells split 2 cell types, 100 genes, gene 1 mostly 0s

  mat_dense <- matrix(sample(1:100, 10000, replace = TRUE), nrow = 100)
  mat_dense[1, 1:95] <- 0
  rownames(mat_dense) <- paste0("Gene", 1:100)
  colnames(mat_dense) <- paste0("Cell", 1:100)
  mat_sparse <- Matrix(mat_dense, sparse = TRUE)


  meta <- data.frame(
    ID = paste0("Cell", 1:100),
    Cell_type = rep(c("Type1", "Type2"), each = 50),
    stringsAsFactors = FALSE
  )

  pc_df <- data.frame(Symbol = paste0("Gene", 1:100))


  return(list(
    mat_dense = mat_dense,
    mat_sparse = mat_sparse,
    meta = meta,
    pc_df = pc_df
  ))

}




# colrank_mat(): focus on default ties/NA arguments
# ------------------------------------------------------------------------------


# Test for basic functionality
test_that("colrank_mat ranks columns correctly", {

  mock_data <- matrix(c(3, 1, 4, 2, 5, 9, 8, 6, 7), nrow = 3, ncol = 3)
  expected_output <- matrix(c(2, 3, 1, 3, 2, 1, 1, 3, 2), nrow = 3, ncol = 3)
  result <- colrank_mat(mock_data)

  expect_equal(result, expected_output)
})



# Test for handling NA values
test_that("colrank_mat handles NA values correctly", {

  mock_data <- matrix(c(3, NA, 4, 2, 5, NA, 8, 6, 7), nrow = 3, ncol = 3)
  expected_output <- matrix(c(2, NA, 1, 2, 1, NA, 1, 3, 2), nrow = 3, ncol = 3)
  result <- colrank_mat(mock_data)

  expect_equal(result, expected_output)
})



# Test for assigning ties
test_that("colrank_mat assigns ties to minimum values", {

  mock_data <- matrix(c(2, 0, 0, -1, 3, 1, 4, 3, 3, 3, 4, 4), nrow = 4, ncol = 3)
  expected_output <- matrix(c(1, 2, 2, 4, 2, 4, 1, 2, 3, 3, 1, 1), nrow = 4, ncol = 3)
  result <- colrank_mat(mock_data)

  expect_equal(result, expected_output)
})



# allrank_mat(): focus on default ties/arg arguments
# ------------------------------------------------------------------------------


# Test for basic functionality
test_that("allrank_mat ranks columns correctly", {

  mock_data <- matrix(c(3, 1, 4, 2, 5, 9, 8, 6, 7), nrow = 3, ncol = 3)
  expected_output <- matrix(c(7, 9, 6, 8, 5, 1, 2, 4, 3), nrow = 3, ncol = 3)
  result <- allrank_mat(mock_data)

  expect_equal(result, expected_output)
})



# Test for handling NA values
test_that("allrank_mat handles NA values correctly", {

  mock_data <- matrix(c(3, NA, 4, 2, 5, NA, 8, 6, 7), nrow = 3, ncol = 3)
  expected_output <- matrix(c(6, NA, 5, 7, 4, NA, 1, 3, 2), nrow = 3, ncol = 3)
  result <- allrank_mat(mock_data)

  expect_equal(result, expected_output)
})



# Test for assigning ties
test_that("allrank_mat assigns ties to minimum values", {

  mock_data <- matrix(c(-1, 0, 6, 3, 0, 4, -2, 0, 4), nrow = 3, ncol = 3)
  expected_output <- matrix(c(8, 5, 1, 4, 5, 2, 9, 5, 2), nrow = 3, ncol = 3)
  result <- allrank_mat(mock_data)

  expect_equal(result, expected_output)
})



# fisherz()
# ------------------------------------------------------------------------------


test_that("fisherz works as expected", {

  mat <- matrix(c(1, 0.54, 0.54, 1), 2, 2)
  result <- fisherz(mat)
  # expected <- DescTools::FisherZ(mat)
  expected <- matrix(c(Inf, 0.6041556, 0.6041556, Inf), 2, 2)

  expect_equal(result, expected, tolerance = 1e-6)
  expect_equal(fisherz(0), 0)
  expect_equal(fisherz(1), Inf)

  expect_error(fisherz(-1.01))
  expect_error(fisherz(1.01))
  expect_error(fisherz(NA))

})



# sparse_pcor()
# ------------------------------------------------------------------------------


test_that("sparse_pcor works correctly with valid input", {

  set.seed(1)

  mat_sparse <- Matrix::rsparsematrix(10, 10, density = 0.5)
  mat_sparse[, 1] <- 0
  colnames(mat_sparse) <- paste0("gene", 1:10)

  mat_dense <- as.matrix(mat_sparse)

  result_sparse <- sparse_pcor(mat_sparse)

  result_dense <-
    suppressWarnings(cor(mat_dense, method = "pearson", use = "pairwise.complete.obs"))

  expect_equal(result_sparse, result_dense, tolerance = 1e-6)

  # All 0 calls should result in NA cors, including self, and not 1
  expect_true(any(is.na(result_sparse)))
  expect_true(is.na(result_sparse[1, 1]))


})




test_that("sparse_pcor throws an error with invalid input", {

  mat_dense <- matrix(c(1, 0, 3, 0, 5, 6, 0, 8, 9), nrow = 3)

  expect_error(sparse_pcor(mat_dense))

})



# calc_sparse_correlation()
# ------------------------------------------------------------------------------


test_that("calc_sparse_correlation throws an error for invalid cor_method", {

  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)

  expect_error(calc_sparse_correlation(mat, cor_method = "invalid"))

})



# Matrix helpers
# ------------------------------------------------------------------------------


test_that("na_to_zero works correctly", {

  mat <- matrix(c(1, NA, 3, 4), 2, 2)
  result <- na_to_zero(mat)
  expected <- matrix(c(1, 0, 3, 4), 2, 2)

  expect_equal(result, expected)

})



test_that("diag_to_one works correctly", {

  mat <- matrix(c(1, 2, 3, 4), 2, 2)
  result <- diag_to_one(mat)
  expected <- matrix(c(1, 2, 3, 1), 2, 2)

  expect_equal(result, expected)

})



test_that("uppertri_to_na works correctly", {

  mat <- matrix(c(1, 2, 3, 4), 2, 2)
  result <- uppertri_to_na(mat)
  expected <- matrix(c(1, 2, NA, 4), 2, 2)

  expect_equal(result, expected)

})



test_that("lowertri_to_symm works correctly", {

  mat <- matrix(c(1, 2, 3, 4), 2, 2)
  result <- lowertri_to_symm(mat)
  expected <- matrix(c(1, 2, 2, 4), 2, 2)

  expect_equal(result, expected)

})



# zero_sparse_cols()
# ------------------------------------------------------------------------------


test_that("zero_sparse_cols sets columns to 0 correctly for sparse matrix", {

  mat <- Matrix(c(1, 0, 3, 0, 5, 0, 7, 0, 0, 0), nrow = 5, ncol = 2, sparse = TRUE)
  result <- zero_sparse_cols(mat, min_cell = 2)

  expect_true(all(result[, 2] == 0))
  expect_equal(mat[, 1], result[, 1])

})



test_that("zero_sparse_cols 0s entire matrix if fewer cells than min_cell", {

  mat <- Matrix(c(1, 2, 3, 4, 5, 1, 2, 3, 4, 0), nrow = 5, ncol = 2, sparse = TRUE)
  result <- zero_sparse_cols(mat, min_cell = 6)

  expect_true(all(result == 0))

})



test_that("zero_sparse_cols preserves mat if min_count is 0", {

  mat <- Matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1, sparse = TRUE)
  result <- zero_sparse_cols(mat, min_cell = 0)
  expect_false(any(result[, 1] == 0))

})



test_that("zero_sparse_cols checks arguments", {

  mat_dense <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
  mat_sparse <- Matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1, sparse = TRUE)

  expect_error(zero_sparse_cols(mat_dense, min_cell = 1))
  expect_error(zero_sparse_cols(mat_sparse, min_cell = -1))

})



# prepare_celltype_mat()
# ------------------------------------------------------------------------------


test_that("prepare_celltype_mat works with sparse matrix and default settings", {

  test_data <- generate_test_data()

  result <- prepare_celltype_mat(mat = test_data$mat_sparse,
                                 meta = test_data$meta,
                                 pc_df = test_data$pc_df,
                                 cell_type = "Type1")

  expect_true(all(result[, 1] == 0))
  expect_equal(nrow(result), 50)
  expect_equal(ncol(result), 100)
  expect_true(all(rownames(result) %in% test_data$meta$ID))

})



test_that("prepare_celltype_mat handles case where no genes meet the min_cell", {

  test_data <- generate_test_data()

  mat_dense <- test_data$mat_dense
  mat_dense[2:100, ] <- mat_dense[1, ]
  mat_sparse <- Matrix(mat_dense, sparse = TRUE)

  result_sparse <- prepare_celltype_mat(mat = mat_sparse,
                                        meta = test_data$meta,
                                        pc_df = test_data$pc_df,
                                        cell_type = "Type1")

  expect_true(all(result_sparse == 0))

  # Also check that matrix is all 0 if fewer cells than min count
  result_sparse <- prepare_celltype_mat(mat = test_data$mat_sparse,
                                        meta = test_data$meta,
                                        pc_df = test_data$pc_df,
                                        cell_type = "Type1",
                                        min_cell = nrow(test_data$mat_sparse) + 1)

  expect_true(all(result_sparse == 0))


})



test_that("prepare_celltype_mat handles invalid cell_type", {

  test_data <- generate_test_data()

  expect_error(prepare_celltype_mat(mat = test_data$mat_sparse,
                                    meta = test_data$meta,
                                    pc_df = test_data$pc_df,
                                    cell_type = "InvalidType"))

})



test_that("prepare_celltype_mat handles incorrect metadata format", {

  test_data <- generate_test_data()
  invalid_meta <- data.frame(ID = paste0("Cell", 1:10), stringsAsFactors = FALSE)

  expect_error(prepare_celltype_mat(mat = test_data$mat_sparse,
                                    meta = invalid_meta,
                                    pc_df = test_data$pc_df,
                                    cell_type = "Type1"))

})



test_that("prepare_celltype_mat handles arguments", {

  test_data <- generate_test_data()

  # Error if dense matrix supplied as argument
  expect_error(prepare_celltype_mat(mat = test_data$mat_dense,
                                    meta = test_data$meta,
                                    pc_df = test_data$pc_df,
                                    cell_type = "Type1"))

  # Error with min cell argument
  expect_error(prepare_celltype_mat(mat = test_data$mat_sparse,
                                    meta = test_data$meta,
                                    pc_df = test_data$pc_df,
                                    cell_type = "Type1",
                                    min_cell = -1))
})



test_that("prepare_celltype_mat checks all cell IDs are in matrix", {

  test_data <- generate_test_data()

  # IDs found that are in meta but not count matrix

  add_mock <- data.frame(ID = paste0("Cell", 101:110),
                         Cell_type = rep("Type1", 10))

  meta <- rbind(test_data$meta, add_mock)

  expect_error(prepare_celltype_mat(mat = test_data$mat_sparse,
                                    meta = meta,
                                    pc_df = test_data$pc_df,
                                    cell_type = "Type1"))

  # No common IDs found at all

  meta$ID <- paste0("Cell", 201:310)

  expect_error(prepare_celltype_mat(mat = test_data$mat_sparse,
                                    meta = meta,
                                    pc_df = test_data$pc_df,
                                    cell_type = "Type1"))

})



test_that("prepare_celltype_mat ensures symbol match between matrix and gene table", {

  # Works when genes in count matrix are a subset of those in gene table
  test_data <- generate_test_data()
  test_data$pc_df <- test_data$pc_df[1:50, , drop = FALSE]

  result <- prepare_celltype_mat(mat = test_data$mat_sparse,
                                 meta = test_data$meta,
                                 pc_df = test_data$pc_df,
                                 cell_type = "Type1")

  expect_equal(ncol(result), 50)
  expect_equal(colnames(result), test_data$pc_df$Symbol)

  # Throws an error when not all genes in gene table are in count matrix
  test_data$pc_df <- rbind(test_data$pc_df,
                           data.frame(Symbol = paste0("Gene", 200:230)))

  expect_error(prepare_celltype_mat(mat = test_data$mat_sparse,
                                    meta = test_data$meta,
                                    pc_df = test_data$pc_df,
                                    cell_type = "Type1"))

})




# init_agg_mat()
# ------------------------------------------------------------------------------


# Initiate a matrix of 0s for holding aggregate correlations
test_that("init_agg_mat works correctly", {

  pc_df <- data.frame(Symbol = c("Gene1", "Gene2", "Gene3"))
  pc_df_dupl <- data.frame(Symbol = c("Gene1", "Gene2", "Gene3", "Gene3"))

  result <- init_agg_mat(pc_df)
  expected <- matrix(0, 3, 3)
  rownames(expected) <- colnames(expected) <- c("Gene1", "Gene2", "Gene3")

  expect_equal(result, expected)

  # Error if pc_df has non-unique entries
  expect_error(init_agg_mat(pc_df_dupl))

})



# increment_na_mat()
# ------------------------------------------------------------------------------


# Count the NAs in cmat and increment the corresponding indices of na_mat
test_that("increment_na_mat works correctly", {

  cmat <- matrix(c(NA, 0.2, 0.3, NA), 2, 2)
  na_mat <- matrix(0, 2, 2)
  result <- increment_na_mat(cmat, na_mat)
  expected <- matrix(c(1, 0, 0, 1), 2, 2)

  expect_equal(result, expected)

})



# transform_correlation_mat()
# ------------------------------------------------------------------------------


test_that("transform_correlation_mat works correctly with allrank", {

  mat <- matrix(c(
    1, 0.8, 0.3,
    0.8, 1, 0.5,
    0.3, 0.5, 1
  ), nrow = 3, ncol = 3)

  rownames(mat) <- colnames(mat) <- c("gene1", "gene2", "gene3")
  result <- transform_correlation_mat(mat, "allrank")

  expected <- matrix(c(
    1, 0.75, 0.25,
    NA, 1, 0.5,
    NA, NA, 1
  ), nrow = 3, ncol = 3)

  rownames(expected) <- colnames(expected) <- c("gene1", "gene2", "gene3")

  expect_equal(result, expected, tolerance = 1e-6)

})



test_that("transform_correlation_mat works correctly with colrank", {

  mat <- matrix(c(
    1, 0.8, -0.1,
    0.8, 1, 0.5,
    -0.1, 0.5, 1
  ), nrow = 3, ncol = 3)

  rownames(mat) <- colnames(mat) <- c("gene1", "gene2", "gene3")
  result <- transform_correlation_mat(mat, "colrank")

  expected <- matrix(c(
    1, (2/3), (1/3),
    (2/3), 1, (1/3),
    (1/3), (2/3), 1
  ), nrow = 3, ncol = 3)

  rownames(expected) <- colnames(expected) <- c("gene1", "gene2", "gene3")

  expect_equal(result, expected, tolerance = 1e-6)

})



test_that("transform_correlation_mat works correctly with FZ", {

  mat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2")
  result <- transform_correlation_mat(mat, "FZ")
  expected <- matrix(c(Inf, 0.5493061, 0.5493061, Inf), 2, 2)
  rownames(expected) <- colnames(expected) <- c("gene1", "gene2")

  expect_equal(result, expected, tolerance = 1e-6)

})



test_that("transform_correlation_mat fails with incorrect arguments", {

  mat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2")
  mat_df <- data.frame(mat)
  mat_badname <- mat
  rownames(mat_badname) <- rev(rownames(mat))

  expect_error(transform_correlation_mat(mat, "invalid_method"))
  expect_error(transform_correlation_mat(mat_df, "FZ"))
  expect_error(transform_correlation_mat(mat_badname, "FZ"))

})



test_that("transform_correlation_mat handles NA values correctly", {

  mat <- matrix(c(
    1, 0.8, NA,
    0.8, 1, 0.5,
    NA, 0.5, 1
  ), nrow = 3, ncol = 3)

  rownames(mat) <- colnames(mat) <- c("gene1", "gene2", "gene3")

  result_allrank <- transform_correlation_mat(mat, "allrank")
  result_colrank <- transform_correlation_mat(mat, "colrank")
  result_fz <- transform_correlation_mat(mat, "FZ")


  expected_allrank <- matrix(c(
    1, 0.75, 0.25,
    NA, 1, 0.5,
    NA, NA, 1
  ), nrow = 3, ncol = 3)

  expected_colrank <- matrix(c(
    1, (2/3), (1/3),
    (2/3), 1, (1/3),
    (1/3), (2/3), 1
  ), nrow = 3, ncol = 3)

  expected_fz <- matrix(c(
    Inf, 1.098612, 0,
    1.098612, Inf, 0.5493061,
    0, 0.5493061, Inf
  ), nrow = 3, ncol = 3)


  rownames(expected_allrank) <- colnames(expected_allrank) <- c("gene1", "gene2", "gene3")
  rownames(expected_colrank) <- colnames(expected_colrank) <- c("gene1", "gene2", "gene3")
  rownames(expected_fz) <- colnames(expected_fz) <- c("gene1", "gene2", "gene3")

  expect_equal(result_allrank, expected_allrank, tolerance = 1e-6)
  expect_equal(result_colrank, expected_colrank, tolerance = 1e-6)
  expect_equal(result_fz, expected_fz, tolerance = 1e-6)

})



# finalize_agg_mat()
# ------------------------------------------------------------------------------


test_that("finalize_agg_mat works correctly with allrank", {

  mat <- matrix(c(
    1, 0.75, 0.25,
    NA, 1, 0.5,
    NA, NA, 1
  ), nrow = 3, ncol = 3) * 2 # as if summing 2 rank std. matrices

  rownames(mat) <- colnames(mat) <- c("gene1", "gene2", "gene3")

  na_mat <- matrix(0, nrow = 3, ncol = 3)
  n_cts <- 2

  result <- finalize_agg_mat(mat, agg_method = "allrank", n_celltypes =  n_cts, na_mat = na_mat)

  expected <- matrix(c(
    1, 0.75, 0.25,
    0.75, 1, 0.5,
    0.25, 0.5, 1
  ), nrow = 3, ncol = 3)

  rownames(expected) <- colnames(expected) <- c("gene1", "gene2", "gene3")

  expect_equal(result, expected, tolerance = 1e-6)

})



test_that("finalize_agg_mat works correctly with colrank", {

  mat <- matrix(c(
    1, (2/3), (1/3),
    (2/3), 1, (1/3),
    (1/3), (2/3), 1
  ), nrow = 3, ncol = 3) * 2 # as if summing 2 rank std. matrices

  rownames(mat) <- colnames(mat) <- c("gene1", "gene2", "gene3")

  na_mat <- matrix(0, nrow = 3, ncol = 3)
  n_cts <- 2

  result <- finalize_agg_mat(mat, agg_method = "colrank", n_celltypes =  n_cts, na_mat = na_mat)

  expected <- matrix(c(
    1, (2/3), (1/3),
    (2/3), 1, (1/3),
    (1/3), (2/3), 1
  ), nrow = 3, ncol = 3)

  rownames(expected) <- colnames(expected) <- c("gene1", "gene2", "gene3")

  expect_equal(result, expected, tolerance = 1e-6)
})



test_that("finalize_agg_mat works correctly with FZ", {

  # Assume we have 3 cell types (n_cts = 3)
  # and the following Fisher's Z transformed values for single measurements:
  # fisherz(0.8) ~ 1.098612, fisherz(0.5) ~ 0.5493061, fisherz(0.6) ~ 0.693147.
  #
  # We simulate an aggregated matrix where:
  # - gene1-gene2 was measured in 2 cell types (na_mat = 1)
  #   => aggregated value = 2 * 1.098612 = 2.197224.
  # - gene1-gene3 was measured in all 3 cell types (na_mat = 0)
  #   => aggregated value = 3 * 0.5493061 = 1.6479183.
  # - gene2-gene3 was measured in only 1 cell type (na_mat = 2)
  #   => aggregated value = 1 * 0.693147 = 0.693147.
  # For the diagonals (self correlations), Fisher's Z of 1 is Inf.

  n_cts <- 3

  mat <- matrix(c(
    Inf,         2.197224,   1.6479183,
    2.197224,    Inf,        0.693147,
    1.6479183,   0.693147,   Inf
  ), nrow = 3, ncol = 3)
  rownames(mat) <- colnames(mat) <- c("gene1", "gene2", "gene3")

  # NA matrix tracking missing measurements:
  # - gene1-gene2: 1 missing (so measured count = 2)
  # - gene1-gene3: 0 missing (measured count = 3)
  # - gene2-gene3: 2 missing (measured count = 1)
  na_mat <- matrix(c(
    0, 1, 0,
    1, 0, 2,
    0, 2, 0
  ), nrow = 3, byrow = TRUE)
  rownames(na_mat) <- colnames(na_mat) <- c("gene1", "gene2", "gene3")

  result <- finalize_agg_mat(mat, agg_method = "FZ", n_celltypes = n_cts, na_mat = na_mat)

  # The expected output divides each element by the count of times it was msrd
  expected <- matrix(c(
    Inf, 2.197224 / (3 - 1), 1.6479183 / (3 - 0),
    2.197224 / (3 - 1), Inf, 0.693147 / (3 - 2),
    1.6479183 / (3 - 0), 0.693147 / (3 - 2), Inf
  ), nrow = 3, ncol = 3)
  rownames(expected) <- colnames(expected) <- c("gene1", "gene2", "gene3")

  expect_equal(result, expected, tolerance = 1e-6)
})



# aggr_coexpr_single_dataset()
# ------------------------------------------------------------------------------



test_that("aggr_coexpr_single_dataset works for Pearson correlation", {

  test_data <- generate_test_data()

  result_allrank <- aggr_coexpr_single_dataset(mat = test_data$mat_sparse,
                                               meta = test_data$meta,
                                               pc_df = test_data$pc_df,
                                               cor_method = "pearson",
                                               agg_method = "allrank",
                                               verbose = FALSE)

  result_colrank <- aggr_coexpr_single_dataset(mat = test_data$mat_sparse,
                                               meta = test_data$meta,
                                               pc_df = test_data$pc_df,
                                               cor_method = "pearson",
                                               agg_method = "colrank",
                                               verbose = FALSE)

  result_fz <- aggr_coexpr_single_dataset(mat = test_data$mat_sparse,
                                          meta = test_data$meta,
                                          pc_df = test_data$pc_df,
                                          cor_method = "pearson",
                                          agg_method = "FZ",
                                          verbose = FALSE)

  expect_true(is.matrix(result_allrank$Agg_mat))
  expect_true(is.matrix(result_colrank$Agg_mat))
  expect_true(is.matrix(result_fz$Agg_mat))

  expect_true(is.matrix(result_allrank$NA_mat))
  expect_equal(dim(result_allrank$Agg_mat), dim(result_allrank$NA_mat))
  expect_equal(result_allrank$NA_mat, result_colrank$NA_mat)
  expect_equal(result_allrank$NA_mat, result_fz$NA_mat)

  expect_equal(dim(result_allrank$Agg_mat),
               c(nrow(test_data$mat_sparse), nrow(test_data$mat_sparse)))

  expect_equal(dim(result_allrank$Agg_mat), dim(result_colrank$Agg_mat))
  expect_equal(dim(result_allrank$Agg_mat), dim(result_fz$Agg_mat))

  expect_equal(diag(result_allrank$Agg_mat), diag(result_colrank$Agg_mat))

})



test_that("aggr_coexpr_single_dataset performs basic argument checks", {

  test_data <- generate_test_data()

  expect_error(aggr_coexpr_single_dataset(mat = NULL,
                                          meta = test_data$meta,
                                          pc_df = test_data$pc_df,
                                          cor_method = "pearson",
                                          agg_method = "allrank",
                                          verbose = FALSE))

  expect_error(aggr_coexpr_single_dataset(mat = test_data$mat_sparse,
                                          meta = NULL,
                                          pc_df = test_data$pc_df,
                                          cor_method = "pearson",
                                          agg_method = "allrank",
                                          verbose = FALSE))

  expect_error(aggr_coexpr_single_dataset(mat = test_data$mat_sparse,
                                          meta = test_data$meta,
                                          pc_df = test_data$pc_df,
                                          # cor_method = "pearson",
                                          agg_method = "allrank",
                                          verbose = FALSE))

  expect_error(aggr_coexpr_single_dataset(mat = test_data$mat_sparse,
                                          meta = test_data$meta,
                                          pc_df = test_data$pc_df,
                                          cor_method = "pearson",
                                          # agg_method = "allrank",
                                          verbose = FALSE))



})



# test_that("aggr_coexpr_single_dataset handles insufficient counts", {
#
#   # Set min_cell to a high number to force insufficient counts
#   result <- aggr_coexpr_single_dataset(mat = test_data$mat_sparse,
#                                        meta = test_data$meta,
#                                        pc_df = test_data$pc_df,
#                                        cor_method = "pearson",
#                                        agg_method = "allrank",
#                                        min_cell = 51)
#
#   expect_true(is.matrix(result$Agg_mat))
#   expect_true(is.matrix(result$NA_mat))
#
#   expect_equal(dim(result$NA_mat), c(10, 10))
#   expect_true(all(result$NA_mat == 1))
#
# })



# load_scdat()
# ------------------------------------------------------------------------------


test_that("load_scdat loads valid paths and checks object", {

  #  Generate mock data and give list expected names before saving temp .RDS
  test_data <- generate_test_data()
  test_data <- list(Mat = test_data$mat_sparse, Meta = test_data$meta)
  temp_file <- tempfile(fileext = ".RDS")
  saveRDS(test_data, temp_file)

  # Mock data with cell IDs that do not match between metadata and matrix
  mismatch_test_data <- test_data
  colnames(mismatch_test_data$Mat) <- as.character(sample(1:100, 100))
  mismatch_temp_file <- tempfile(fileext = ".RDS")
  saveRDS(mismatch_test_data, mismatch_temp_file)

  result <- load_scdat(temp_file)

  expect_true(inherits(result$Mat, "dgCMatrix"))
  expect_true(inherits(result$Meta, "data.frame"))

  expect_equal(result$Mat, test_data$Mat)
  expect_equal(result$Meta, test_data$Meta)

  expect_error(load_scdat(path = "A/surely/fake/path.RDS"))
  expect_error(load_scdat(mismatch_temp_file))

  unlink(temp_file)
  unlink(mismatch_temp_file)

})




# aggr_coexpr_multi_dataset()
# ------------------------------------------------------------------------------



test_that("aggr_coexpr_multi_dataset works for Pearson correlation", {

  #  Generate mock data to load inside call
  test_data <- generate_test_data()
  test_data1 <- list(Mat = test_data$mat_sparse, Meta = test_data$meta)
  test_data2 <- test_data1
  temp_file1 <- tempfile(fileext = ".RDS")
  temp_file2 <- tempfile(fileext = ".RDS")
  saveRDS(test_data1, temp_file1)
  saveRDS(test_data2, temp_file2)

  input_df <- data.frame(ID = c("ID1", "ID2"),
                         Path = c(temp_file1, temp_file2),
                         Cell_type = c("Type1", "Type1"))


  result_allrank <- aggr_coexpr_multi_dataset(input_df = input_df,
                                              pc_df = test_data$pc_df,
                                              cor_method = "pearson",
                                              agg_method = "allrank",
                                              verbose = FALSE)

  result_colrank <- aggr_coexpr_multi_dataset(input_df = input_df,
                                              pc_df = test_data$pc_df,
                                              cor_method = "pearson",
                                              agg_method = "colrank",
                                              verbose = FALSE)

  result_fz <- aggr_coexpr_multi_dataset(input_df = input_df,
                                         pc_df = test_data$pc_df,
                                         cor_method = "pearson",
                                         agg_method = "FZ",
                                         verbose = FALSE)

  expect_true(is.matrix(result_allrank$Agg_mat))
  expect_true(is.matrix(result_colrank$Agg_mat))
  expect_true(is.matrix(result_fz$Agg_mat))

  expect_true(is.matrix(result_allrank$NA_mat))
  expect_equal(dim(result_allrank$Agg_mat), dim(result_allrank$NA_mat))
  expect_equal(result_allrank$NA_mat, result_colrank$NA_mat)
  expect_equal(result_allrank$NA_mat, result_fz$NA_mat)

  expect_equal(dim(result_allrank$Agg_mat),
               c(nrow(test_data$mat_sparse), nrow(test_data$mat_sparse)))

  expect_equal(dim(result_allrank$Agg_mat), dim(result_colrank$Agg_mat))
  expect_equal(dim(result_allrank$Agg_mat), dim(result_fz$Agg_mat))

  expect_equal(diag(result_allrank$Agg_mat), diag(result_colrank$Agg_mat))


  unlink(temp_file1)
  unlink(temp_file2)

})

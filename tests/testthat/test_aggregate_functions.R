# TODO: more explicit check of NA content in sparse_cor()


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



# zero_sparse_cols()
# ------------------------------------------------------------------------------


test_that("zero_sparse_cols sets columns to 0 correctly for sparse matrix", {

  mat <- Matrix(c(1, 0, 3, 0, 5, 0, 7, 0, 0, 0), nrow = 5, ncol = 2, sparse = TRUE)
  result <- zero_sparse_cols(mat, min_count = 2)

  expect_true(all(result[, 2] == 0))
  expect_equal(mat[, 1], result[, 1])

})



test_that("zero_sparse_cols handles edge cases for sparse matrix", {

  mat <- Matrix(c(1, 2, 3, 4, 5, 1, 2, 3, 4, 0), nrow = 5, ncol = 2, sparse = TRUE)
  result <- zero_sparse_cols(mat, min_count = 5)

  expect_true(all(result[, 2] == 0))
  expect_error(zero_sparse_cols(mat, min_count = 6))

})



test_that("zero_sparse_cols handles minimum count of 0 for sparse matrix", {

  mat <- Matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1, sparse = TRUE)
  result <- zero_sparse_cols(mat, min_count = 0)
  expect_false(any(result[, 1] == 0))

})



test_that("zero_sparse_cols checks arguments", {

  mat_dense <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
  mat_sparse <- Matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1, sparse = TRUE)

  expect_error(zero_sparse_cols(mat_dense, min_count = 1))
  expect_error(zero_sparse_cols(mat_sparse, min_count = -1))
  expect_error(zero_sparse_cols(mat_sparse, min_count = 10))

})

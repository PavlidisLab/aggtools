# TODO: more explicit check of NA content in sparse_pcor()


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

})




test_that("sparse_pcor throws an error with invalid input", {

  mat_dense <- matrix(c(1, 0, 3, 0, 5, 6, 0, 8, 9), nrow = 3)

  expect_error(sparse_pcor(mat_dense))

})


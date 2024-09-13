# Test for basic functionality
test_that("colrank_mat() ranks columns correctly", {

  mock_data <- matrix(c(3, 1, 4, 2, 5, 9, 8, 6, 7), nrow = 3, ncol = 3)
  expected_output <- matrix(c(2, 3, 1, 3, 2, 1, 1, 3, 2), nrow = 3, ncol = 3)
  result <- colrank_mat(mock_data)

  expect_equal(result, expected_output)
})



# Test for handling NA values
test_that("colrank_mat() handles NA values correctly", {

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

# Load required libraries
library(testthat)

# Test for gaussian change in mean where pre-change mean is known
test_that("FOCuS detects change in mean when pre-change mean is known", {
  set.seed(123)
  # Generate data with a change in mean at position 50
  data <- c(rnorm(5000, mean = 0), rnorm(50, mean = 1))

  # Run FOCuS on the data
  result <- FOCuS(data, thres = 15, mu0 = 0)

  # Check that the detected changepoint is close to the true changepoint
  expect_true(abs(result$changepoint - 5000) <= 5)
})

# Test for gaussian change in mean where pre-change mean is unknown
test_that("FOCuS detects change in mean when pre-change mean is unknown", {
  set.seed(123)
  # Generate data with a change in mean at position 50
  data <- c(rnorm(5000, mean = 0), rnorm(50, mean = 1))

  # Run FOCuS on the data
  result <- FOCuS(data, thres = 15)

  # Check that the detected changepoint is close to the true changepoint
  expect_true(abs(result$changepoint - 5000) <= 5)
})

# Test for gaussian change in mean where pre-change mean is known but there are outliers
test_that("FOCuS detects change in mean when pre-change mean is known and there are outliers", {
  set.seed(123)
  # Generate data with a change in mean at position 50 and add some outliers
  data <- c(rnorm(5000, mean = 0), rnorm(50, mean = 1))
  data[c(10,20)] <- c(10,-10)

  # Run FOCuS on the data with K set to ignore outliers larger than K
  result <- FOCuS(data, thres = 15, mu0 = 0, K = 3)

  # Check that the detected changepoint is close to the true changepoint
  expect_true(abs(result$changepoint - 5000) <=5)
})
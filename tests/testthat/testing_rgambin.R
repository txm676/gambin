test_that("Testing dgambin", {
  skip_on_cran()
  
  N = 1e3
  ## Test values outside 0 and maxoctave
  x = rgambin(N, 1, 1, w = 20)
  testthat::expect_length(x, N)
  testthat::expect_true(max(x) <= 1)
  testthat::expect_true(min(x) >= 0)
  
  # Test weightings
  x = rgambin(N, c(1,10), c(1,10),  w = c(1, 0))
  testthat::expect_length(x, N)
  testthat::expect_true(max(x) <= 1)
  testthat::expect_true(min(x) >= 0)
  
  # Test weightings
  x = rgambin(N, c(1,10), c(1,10),  w = c(0, 1))
  testthat::expect_length(x, N)
  testthat::expect_true(median(x) >= 3)
  testthat::expect_true(max(x) <= 10)
}
)
test_that("Testing dgambin", {
  skip_on_cran()
  x = dgambin(-1:2, 1, 1)
  testthat::expect_length(x, 4)
  testthat::expect_true((x[1] == 0) && (x[4] == 0))
  
}
)
test_that("Testing dgambin", {
  skip_on_cran()
  
  ## Test values outside 0 and maxoctave
  x = dgambin(-1:2, 1, 1)
  testthat::expect_length(x, 4)
  testthat::expect_true((x[1] == 0) && (x[4] == 0))
  
  testthat::expect_true(x[3] == dgambin(1, 1, 1) && x[4] == dgambin(2, 1, 1))

  ## Test weights
  expect_equivalent(dgambin(0:3, 1, w = 1, maxoctave = 3),
                    dgambin(0:3, 1, maxoctave = 3) )
    
  ## Equal weightings
  expect_equivalent(dgambin(0:3, c(1,1), w = c(1, 1), maxoctave = c(3, 3)),
                    dgambin(0:3, 1, maxoctave = 3) )
  

  # zero weighting
  expect_equivalent(dgambin(0:3, c(1, 2), w = c(0.5, 0), maxoctave = c(3, 4)),
                    dgambin(0:3, 1, maxoctave = 3) )
  
  
  ## Test qgambin
  expect_equivalent(qgambin(0.5, 2, maxoctave = 10),
                    qgambin(0.5, c(2, 2), w=c(1, 1), maxoctave = c(10, 10)))

  
  }
)
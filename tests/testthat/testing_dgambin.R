test_that("Testing dgambin", {
  skip_on_cran()
  x = dgambin(-1:2, 1, 1)
  testthat::expect_length(x, 4)
  testthat::expect_true((x[1] == 0) && (x[4] == 0))

  ## Test weigths
  expect_equivalent(dgambin(0:3, 1, w = 1, maxoctave = 3),
                    dgambin(0:3, 1, maxoctave = 3) )
    
  expect_equivalent(dgambin(0:3, c(1,1), w = c(1, 1), maxoctave = 3),
                    dgambin(0:3, 1, maxoctave = 3) )
  
  expect_equivalent(dgambin(0:3, c(1,1), w = c(0.5, 0.5), maxoctave = 3),
                    dgambin(0:3, 1, maxoctave = 3) )
  
  ## Test qgambin
  expect_equivalent(qgambin(0.5, 2, maxoctave = 10),
                    qgambin(0.5, c(2, 2), w=c(1, 1), maxoctave = 10))

  
  }
)
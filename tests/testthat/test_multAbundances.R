context("mult_abundances")

test_that("fit_abudances returns correct info", {
  skip_on_cran()
  
  mult <- matrix(0, nrow = 50, ncol = 10)
  mult <- apply(mult, 2, function(x) ceiling(rlnorm(length(x), 0, 2.5)))
  mm <- suppressWarnings(mult_abundances(mult, N = 10, subsample = NULL))
  
  expect_equal(length(mm), 4)
  expect_equal(length(mm$Mean.Stan.Alpha), 10)
  expect_is(mm, c("list", "vector"))
  expect_equal(dim(mm$All.Alpha), c(10, 10))
  expect_error(mult_abundances(c(4:6)), "Data are not in correct format")
  
})



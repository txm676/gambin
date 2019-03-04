context("fit_abundances")

test_that("fit_abundances returns correct info", {
  skip_on_cran()
  data(moths)
  f <- fit_abundances(moths)
  f2 <- fit_abundances(moths, no_of_components = 2)
  f3 <- fit_abundances(moths, subsample = 200)
  
  expect_message(fit_abundances(moths, no_of_components = 2), "Using 1 core. Your machine has 12 available.")
  
  expect_equal(round(f$alpha, 3), 1.645)
  expect_equal(round(f$octaves, 0), 10)
  expect_is(f, "gambin")
  
  expect_equal(round(f2$alpha, 2), c(1.92, 5.04))
  expect_equal(round(f2$w, 2), c(0.76, 0.24))
  expect_is(f2, "gambin")
  
  expect_match(f3$Dataname, "200 individuals sampled from moths")
  
  expect_error(fit_abundances(moths, subsample = 50000), "cannot take a sample larger than the population when 'replace = FALSE'")
 
})


test_that("summary_fit_abundances works properly", {
  skip_on_cran()
  data(moths)
  f <- summary(fit_abundances(moths))
  f2 <- summary(fit_abundances(moths, no_of_components = 2))

  expect_is(f, "summary.gambin")
  expect_equal(round(f$ConfInt95, 2), c(1.32, 2.03))
  expect_equal(round(f2$ChiSq$statistic, 2), c("X-squared" = 2.92))
  
  
})
context("deconstruct_modes")

test_that("deconstruct_modes returns correct info", {
  skip_on_cran()
  data(categ)
  fits2 = fit_abundances(categ$abundances, no_of_components = 2)
  d3 <- deconstruct_modes(fits2, dat = categ, peak_val = NULL, abundances = "abundances", 
                          species = "species", categ = NULL, plot_modes = FALSE)
  expect_equal(d3$Peak_locations, c(0, 5))
  
  d4 <- deconstruct_modes(fits2, dat = categ, peak_val = NULL, abundances = "abundances", 
                          species = "species", categ = "status", plot_modes = FALSE)
  expect_equal(d4$Peak_locations, c(0, 5))
  expect_equal(d4$Summary_table[1,2], 74)
  data(moths)
  expect_error(deconstruct_modes(moths), "fit is not a gambin object")
  
})

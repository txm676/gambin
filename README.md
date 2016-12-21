
<!-- README.md is generated from README.Rmd. Please edit that file -->
The `gambin` R package
======================

[![Build Status](https://travis-ci.org/mkborregaard/gambin.svg?branch=master)](https://travis-ci.org/mkborregaard/gambin) [![Downloads](http://cranlogs.r-pkg.org/badges/gambin?color=brightgreen)](https://cran.r-project.org/package=gambin) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/gambin)](https://cran.r-project.org/package=gambin) [![codecov.io](https://codecov.io/github/mkborregaard/gambin/coverage.svg?branch=master)](https://codecov.io/github/mkborregaard/gambin?branch=master)

This package fits the `gambin` distribution to species-abundance distributions from ecological data. 'gambin' is short for 'gamma-binomial'. The main function is `fitGambin`, which estimates the `alpha` parameter of the gambin distribution using maximum likelihood. Functions are also provided to generate the gambin distribution and for calculating likelihood statistics.

The package is currently on CRAN and can be installed via

``` r
install.packages("gambin")
```

Alternatively, you can install the developmental version via

``` r
devtools::install_github("mkborregaard/gambin")
```

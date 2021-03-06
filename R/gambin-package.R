#' @title Fit the gambin model to species abundance distributions
#' @name gambin-package
#' @aliases gambin-package gambin
#' @docType package
#' @description This package provides functions for fitting unimodal and
#'   multimodal gambin distributions to species-abundance distributions from
#'   ecological data. The main function is \code{fit_abundances()}, which
#'   estimates the 'alpha' parameter(s) of the gambin distribution using maximum
#'   likelihood.
#' @details The gambin distribution is a sample distribution based on a
#'   stochastic model of species abundances, and has been demonstrated to fit
#'   empirical data better than the most commonly used species-abundance models
#'   (see references). Gambin is a stochastic model which combines the gamma
#'   distribution with a binomial sampling method. To fit the gambin
#'   distribution, the abundance data is first binned into octaves. The expected
#'   abundance octave of a species is given by the number of successful
#'   consecutive Bernoulli trials with a given parameter \code{p}. The parameter
#'   \code{p} of species is assumed to distributed according to a gamma
#'   distribution. This approach can be viewed as linking the gamma distribution
#'   with the probability of success in a binomial process with x trials. Use
#'   the \code{fit_abundances()} function to fit the gambin model to a vector of
#'   species abundances, optionally using a subsample of the individuals. Use
#'   the \code{mult_abundances()} function to fit the gambin model to multiple
#'   sites / samples and return the alpha values for each model fit (both the
#'   raw values and the alpha values standardised by the number of
#'   individuals).The package estimates the alpha (shape) parameter with
#'   associated confidence intervals. Methods are provided for plotting the
#'   results, and for calculating the likelihood of fits.
#'
#'   The package now provides functionality to fit multimodal gambin
#'   distributions (i.e. a gambin distribution with more than one mode), and to
#'   deconstruct and examine a multimodal gambin model fit
#'   (\code{deconstruct_modes}).
#'   
#' @references Matthews, T.J., Borregaard, M.K., Ugland, K.I., Borges, P.A.V,
#'   Rigal, F., Cardoso, P. and Whittaker, R.J. (2014) The gambin model provides
#'   a superior fit to species abundance distributions with a single free
#'   parameter: evidence, implementation and interpretation. Ecography 37:
#'   1002-1011.
#'
#'   Matthews, T.J., Borregaard, M.K., Gillespie, C.S., Rigal, F., Ugland, K.I.,
#'   Krüger, R.F., Marques, R., Sadler, J.P., Borges, P.A.V., Kubota, Y. &
#'   Whittaker, R.J. (2019) Extension of the gambin model to multimodal species
#'   abundance distributions. Methods in Ecology and Evolution, 10, 432-437.
#'   
#'   Ugland, K.I., Lambshead, F.J.D., McGill, B.J., Gray, J.S., O'Dea, N.,
#'   Ladle, R.J. & Whittaker, R.J. (2007). Modelling dimensionality in species
#'   abundance distributions: description and evaluation of the Gambin model.
#'   Evolutionary Ecology Research, 9, 313-324.
#' @keywords package
#' @seealso \url{https://github.com/txm676/gambin}
#' @examples 
#' data(moths, package = "gambin")
#' fit = fit_abundances(moths)
#' barplot(fit)
#' lines(fit)
#' AIC(fit)
NULL


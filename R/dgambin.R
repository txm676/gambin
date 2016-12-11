#' @title Calculate the gambin distribution 
#' @description Calculate the expected number of species in octaves for a given value of alpha and maxoctave
#' @param x vector of (non-negative integer) quantiles.
#' @param alpha The shape parameter of the GamBin distribution.
#' @param maxoctave The scale parameter of the GamBin distribution - which octave is the highest in the empirical dataset?
#' @param log logical; Tf \code{TRUE}, probabilities p are given as log(p).
#' @param total_species The total number of species in the empirical dataset
#' @details   \code{dgambin} gives the distribution function of gambin, so all octaves sum to 1.
#' \code{gambin_exp} multiplies this by the total number of species to give the expected GamBin distribution in units of species, 
#' for comparison with empirical data.  
#' @return A vector with length MaxOctave + 1 of the expected number of species in each octave
#' @references  Ugland, K.I., Lambshead, F.J.D., McGill, B.J., Gray, J.S., O'Dea, N., Ladle, R.J. & Whittaker, R.J. (2007). Modelling dimensionality in species abundance distributions: description and evaluation of the Gambin model. Evolutionary Ecology Research, 9, 313-324.
#' Matthews, T.J., Borregaard, M.K., Ugland, K., Borges, P.A.V., Rigal, F. & Whittaker, R.J. (Early View Online). 
#' The gambin model provides a superior fit to species abundance distributions with a single free parameter: evidence, 
#' implementation and interpretation. Ecography, DOI: 10.1111/ecog.00861
#' @examples 
#' expected <- gambin_exp(4, 13, 200)
#' plot(expected, type = "l")
#' @export
dgambin = function(x, alpha, maxoctave, log = FALSE)
{
  vec = 0:100/100
  # Calculates the 'fitness' distribution of species for a given alpha
  qG99 = qgamma(0.99, alpha, 1) * vec
  Gj = (pgamma(qG99[-1], alpha, 1) - pgamma(qG99[-101], alpha, 1)) / 0.99
  
  gambin_p = function(k) {
    if(k < 0 || k > maxoctave) 0
    else sum(choose(maxoctave, k) * vec[-1]^k * (1 - vec[-1])^(maxoctave - k) * Gj)
  }
  
  # Apply Pk to each octave:
  res = vapply(x, gambin_p,  FUN.VALUE = numeric(1))
  if(log)
    res = log(res)
  
  res
}

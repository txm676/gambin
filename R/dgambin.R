#' @title Calculate the gambin distribution 
#' @description Calculate the expected number of species in octaves for a given value of alpha and maxoctave
#' @param alpha The shape parameter of the GamBin distribution
#' @param maxoctave The scale parameter of the GamBin distribution - which octave is the highest in the empirical dataset?
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
dgambin <-
function(alpha, maxoctave)
{
  # calculates the 'fitness' distribution of species for a given alpha
  qG99 = qgamma(0.99, alpha, 1) / 100
  b1 = 0:99 * qG99
  b2 = 1:100 * qG99
  Gj <- (pgamma(b2, alpha, 1) - pgamma(b1, alpha, 1)) / 0.99
  
  # a function to give the probability distribution of the GamBin in 
  Pk <- function(k, octmax) # k is in 0:(nOct-1)
    return(sum(choose(octmax, k) * (1:100/100)^k * (1 - 1:100/100)^(octmax - k) * Gj)) 
  
  # applies Pk to each octave:
  ret <- sapply(0:maxoctave, Pk, octmax = maxoctave)
  names(ret) <- 0:maxoctave
  ret
}

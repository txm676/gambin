dgambin_single = function(x, alpha, maxoctave, log = FALSE)
{
  vec = 0:100/100
  # Calculates the 'fitness' distribution of species for a given alpha
  qG99 = qgamma(0.99, alpha, 1) * vec
  Gj = diff(pgamma(qG99, alpha, 1)) / 0.99
  
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

#' @title The  mixture gambin distribution 
#' @description Calculate the expected number of species in octaves for a given value of alpha and maxoctave
#' @param x vector of (non-negative integer) quantiles.
#' @param alpha The shape parameter of the GamBin distribution.
#' @param maxoctave The scale parameter of the GamBin distribution - which octave is the highest in the empirical dataset?
#' @param log logical; If \code{TRUE}, probabilities p are given as log(p).
#' @param total_species The total number of species in the empirical dataset
#' @param w A vector of weights. Default, a single weight. This vector must of the same length as alpha
#' @details \code{dgambin} gives the distribution function of a mixture gambin, so all octaves sum to 1.
#' \code{gambin_exp} multiplies this by the total number of species to give the expected GamBin distribution in units of species, 
#' for comparison with empirical data.  
#' @return A vector with length MaxOctave + 1 of the expected number of species in each octave
#' @references  Ugland, K.I., Lambshead, F.J.D., McGill, B.J., Gray, J.S., O'Dea, N., Ladle, R.J. & Whittaker, R.J. (2007). Modelling dimensionality in species abundance distributions: description and evaluation of the Gambin model. Evolutionary Ecology Research, 9, 313-324.
#' Matthews, T.J., Borregaard, M.K., Ugland, K., Borges, P.A.V., Rigal, F. & Whittaker, R.J. (Early View Online). 
#' The gambin model provides a superior fit to species abundance distributions with a single free parameter: evidence, 
#' implementation and interpretation. Ecography, DOI: 10.1111/ecog.00861
#' @examples 
#' expected = gambin_exp(4, 13, total_species = 200)
#' plot(expected, type = "l")
#' @export
dgambin = function(x, alpha, maxoctave, w = 1,log = FALSE)
{
  if(any(w < 0)) stop("w must be non-negative", call. = FALSE)
  
  w = if(length(w) == 1) rep(1/length(alpha), length(alpha)) else w = w/sum(w)
  
  res = vapply(seq_along(alpha), 
               function(i) w[i]*dgambin_single(x, alpha[i], maxoctave, log=FALSE), 
               FUN.VALUE = numeric(length(x)))
  
  res = rowSums(res)
  
  if(log)
    res = log(res)
  res
}

#' @param q	vector of quantiles.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @param log.p	logical; if \code{TRUE}, probabilities p are given as log(p).
#' @rdname dgambin
#' @export
pgambin = function(q, alpha, maxoctave, w = 1, lower.tail = TRUE, log.p = FALSE) 
{
  # Form CMF from mass function. Then manipulate as necessary
  probs = dgambin(0:maxoctave, alpha,  maxoctave, w=w)
  
  if(!lower.tail)
    probs = 1- probs
  
  cum_probs = cumsum(probs)
  if(log.p)
    cum_probs = log(cum_probs)
  
  cum_probs[floor(q) + 1]
}

#' @param n	number of random values to return.
#' @rdname dgambin
#' @export
rgambin = function(n, alpha, maxoctave, w=1) 
{
  # Initialise parameters
  if(length(n) > 1L) n = length(n)
  # Form look-up table
  probs = dgambin(0:maxoctave, alpha, maxoctave, w=w)
  
  sample(0:maxoctave, prob = probs, replace=TRUE, size = n)
}

#' @param p vector of probabilities.
#' @rdname dgambin
#' @export
qgambin = function(p, alpha, maxoctave, w=1, lower.tail = TRUE, log.p = FALSE)
{
  # Form CMF from mass function. Then manipulate as necessary
  probs = dgambin(0:maxoctave, alpha, maxoctave, w=w)
  
  # Add on 0 for cut
  probs = c(0, probs)
  
  if(!lower.tail)
    probs = 1 - probs
  
  if(log.p) probs = log(probs)
  cum_probs = cumsum(probs)
  
  ## Use cut and exploit factor
  as.numeric(cut(p, cum_probs)) - 1
}






























sample_abundances  = function(abundances, individuals)
{
  if(!is.numeric(abundances)) stop("abundances must be numeric")
  species <- as.factor(1:length(abundances))
  samplevector <- rep(species, abundances)
  samp <- sample(samplevector, individuals, replace = F)
  return(table(samp))
}



#' @title Reclassify a vector of species' abundances into abundance octaves
#' @description Creates abundance octaves by a log2 transform that doubles the number 
#' of abundance classes within each octave (method 3 of Gray, Bjoergesaeter & Ugland 2006). 
#' Octave 0 contains the number of species with 1 individual, octave 1 the number of species 
#' with 2 or 3 individuals, octave 2 the number of species with 4 to 7 individuals, and so forth.  
#' @param abundances A  numerical vector of species abundances in a community.
#' @param subsample If > 0, the community is subsampled by this number of individuals before creating octaves. This is useful for analyses where \code{alpha} is estimated from a standardized number of individuals.
#' @return A data.frame with two variables: \code{octave} with the name of each octave and \code{species} with the number of species in that octave
#' @references Gray, J.S., Bjoergesaeter, A. & Ugland, K.I. (2006) On plotting species abundance distributions. Journal of Animal Ecology, 75, 752-756.
#' @examples 
#' data(moths)
#' create_octaves(moths)
#' @export
create_octaves = function(abundances, subsample = 0)  
{
  if(subsample > 0) abundances <- sample_abundances(abundances, subsample)
  stopifnot(is.numeric(abundances))
  abundances = abundances[abundances > 0]     # remove zeros
  octs = floor(vapply(abundances, log2, FUN.VALUE = numeric(1)))
  octs = factor(octs, levels = 0:max(octs))   # ensure that all octaves are tabled, even if no species fall in that octave
  ret <- data.frame(table(octs))
  names(ret) = c("octave", "species")
  ret$octave = 0:(nrow(ret)-1)                # octaves are numbered from 0, 1, 2... etc.
  ret
}

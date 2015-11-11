\name{The gambin distribution}
\alias{dgambin}
\alias{gambin_exp}

\title{Calculate the gambin distribution }
\description{ calculate the expected number of species in octaves for a given value of alpha and maxoctave}
\usage{
dgambin(alpha, maxoctave)
gambin_exp(alpha, maxoctave, total_species)
}

\arguments{
  \item{alpha}{ the shape parameter of the GamBin distribution }
  \item{maxoctave}{ the scale parameter of the GamBin distribution - which octave is the highest in the empirical dataset? }
  \item{total_species}{the total number of species in the empirical dataset}
  }
\details{
\code{dgambin} gives the distribution function of gambin, so all octaves sum to 1.
\code{gambin_exp} multiplies this by the total number of species to give the expected GamBin distribution in units of species, for comparison with empirical data.  
}
\value{A vector with length MaxOctave + 1 of the expected number of species in each octave
}
\references{
Ugland, K.I., Lambshead, F.J.D., McGill, B.J., Gray, J.S., O'Dea, N., Ladle, R.J. & Whittaker, R.J. (2007). Modelling dimensionality in species abundance distributions: description and evaluation of the Gambin model. Evolutionary Ecology Research, 9, 313-324.\cr
Matthews, T.J., Borregaard, M.K., Ugland, K., Borges, P.A.V., Rigal, F. & Whittaker, R.J. (Early View Online). The gambin model provides a superior fit to species abundance distributions with a single free parameter: evidence, implementation and interpretation. Ecography, DOI: 10.1111/ecog.00861 }

\author{Michael K. Borregaard, Gabriel Areto, Thomas J. Matthews & Karl I. Ugland}

\examples{
expected <- gambin_exp(4, 13, 200)
plot(expected, type = "l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{gambin-package}
\alias{gambin-package}
\alias{gambin}
\docType{package}
\title{
Fit the gambin model to species abundance distributions
}
\description{
This package provides functions for fitting the gambin distribution to species-abundance distributions from ecological data. The main function is fitGambin, which estimates the 'alpha' parameter of the gambin distribution using maximum likelihood. 
}
\details{
\tabular{ll}{
Package: \tab gambin\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2015-11-11\cr
License: \tab GPL-2\cr
}
The gambin distribution is a sample distribution based on a stochastic model of species abundances, and has been demonstrated to fit empirical data better than the most commonly used species-abundance models (see references). Gambin is a stochastic model which combines the gamma distribution with a binomial sampling method. To fit the gambin distribution, the abundance data is first binned into octaves. The expected abundance octave of a species is given by the number of successfull consecutive Bernoulli trials with a given parameter \code{p}. The parameter \code{p} of species is assumed to distributed according to a gamma distribution. This approach can be viewed as linking the gamma distribution with the probability of success in a binomial process with x trials.
Use the \code{fitGambin} function fit the gambin model to a vector of species abundances, optionally using a subsample of the individuals. The package estimates the alpha (shape) parameter with associated confidence intervals. Methods are provided for plotting the results, and for calculating the likelihood of fits.
}
\author{Michael Krabbe Borregaard, Thomas J. Matthews & Karl I. Ugland\cr
Maintainer: Michael Krabbe Borregaard <mkborregaard@snm.ku.dk>}
\references{
Matthews, T.J., Borregaard, M.K., Ugland, K., Borges, P.A.V., Rigal, F. & Whittaker, R.J. (Early View Online). The gambin model provides a superior fit to species abundance distributions with a single free parameter: evidence, implementation and interpretation. Ecography, DOI: 10.1111/ecog.00861 \cr
Ugland, K.I., Lambshead, F.J.D., McGill, B.J., Gray, J.S., O'Dea, N., Ladle, R.J. & Whittaker, R.J. (2007). Modelling dimensionality in species abundance distributions: description and evaluation of the Gambin model. Evolutionary Ecology Research, 9, 313-324.
}

\keyword{ package }
\examples{
data(moths)
fit <- fitGambin(moths)
plot(fit)
AIC(fit)
}

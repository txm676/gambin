\name{logLik.gambin}
\alias{logLik.gambin}
\alias{AIC.gambin}
\alias{AICc}
\alias{AICc.gambin}
\alias{BIC.gambin}

\title{
Likelihood statistics for the GamBin model
}
\description{
Uses likelihood and information theoretical approaches to reveal the degree of fit of the GamBin model to empirical species abundance distributions.
}
\usage{
\method{logLik}{gambin}(object, ...)
\method{AIC}{gambin}(object, ...)
\method{AICc}{gambin}(object, ...)
\method{BIC}{gambin}(object, ...)
}

\arguments{
  \item{object}{ an object of type \code{gambin} }
  \item{\dots}{ further arguments to pass to the function }
}

\value{logLik returns an R object of type \code{logLik}. The other function return the numerical value of the statistic
}
\references{
Akaike, Hirotugu. "A new look at the statistical model identification." Automatic Control, IEEE Transactions on 19.6 (1974): 716-723.
}
\author{
Michael K. Borregaard & Thomas J. Matthews
}

\seealso{
\code{\link{fitGambin}}
}
\examples{

data(moths)
fit <- fitGambin(moths)
AIC(fit)
}


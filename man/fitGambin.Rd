\name{fitGambin}
\alias{fitGambin}
\alias{print.gambin}
\alias{summary.gambin}
\alias{print.summary.gambin}
\alias{predict.gambin}
\alias{confint.gambin}
\alias{plot.gambin}
\alias{nobs.gambin}


\title{ Fit the GamBin model to a species abundance distribution }
\description{ Uses maximum likelihood methods to fit the GamBin model to binned species abundances. To control for the effect of sample size, the abundances may be subsampled prior to fitting. }
\usage{ 
	fitGambin(abundances, subsample = 0) 
	\method{plot}{gambin}(x, barcol = "grey", barwidth = 1, cex.dots = 1, dotpch = 16, 
	dotcol = par("fg"), line = FALSE, lwd = 1, linecol = par("fg"), ...)
	\method{predict}{gambin}(object, ...)
	\method{confint}{gambin}(object, parm = "alpha", level = 0.95, ...)
}



\arguments{
  \item{abundances}{ Either a vector of abundances of all species in the sample/community; or the result of \code{create_octaves} }
  \item{subsample}{ The number of individuals to sample from the community before fitting the GamBin model. If subsample == 0 the entire community is used }
  \item{object}{ a code{gambin} object created by \code{fitGambin} }
  \item{x}{ a code{gambin} object created by \code{fitGambin}}
  \item{parm}{the parameter to calculate confidence intervals from. Only alpha is implemented}
  \item{level}{ the significance level of the confidence intervals }
  \item{barcol}{ the colour of the bars illustrating the empirical abundance of species in octaves }
  \item{barwidth}{ the width of the bars illustrating the empirical abundance of species in octaves }
  \item{cex.dots}{ the size of the dots illustrating the fitted abundance of species in octaves }
  \item{dotpch}{ the point character of the dots illustrating the fitted abundance of species in octaves }
  \item{dotcol}{ the colour of the dots illustrating the fitted abundance of species in octaves }
  \item{line}{ should the dots be connected with a line? }
  \item{lwd}{ the width of the line connecting dots }
  \item{linecol}{ the colour of the line connecting dots }
  \item{\dots}{ further arguments to pass to \code{barplot}}

}

\details{ The gambin distribution is fit to the number of species in abundance octaves, as specified by the \code{create_octaves} function. Because the shape of species abundance distributions depend on sample size, abundances of different communities should be compared on equally large samples. The sample size can be set by the \code{subsample} parameter. To estimate \code{alpha} from a standardised sample, the function must be run several times; see the examples. 
	The \code{plot} method creates a barplot showing the observed number of species in octaves, with the fitted GamBin distribution shown as black dots.
	}
\value{
 The \code{fitGambin} function returns an object of class \code{gambin}, with the \code{alpha} and \code{MaxOctave} parameters of the GamBin distribution, the likelihood of the fit, and the empirical distribution over octaves.
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Ugland, K.I., Lambshead, F.J.D., McGill, B.J., Gray, J.S., O'Dea, N., Ladle, R.J. & Whittaker, R.J. (2007). Modelling dimensionality in species abundance distributions: description and evaluation of the Gambin model. Evolutionary Ecology Research, 9, 313-324.\cr
Matthews, T.J., Borregaard, M.K., Ugland, K., Borges, P.A.V., Rigal, F. & Whittaker, R.J. (Early View Online). The gambin model provides a superior fit to species abundance distributions with a single free parameter: evidence, implementation and interpretation. Ecography, DOI: 10.1111/ecog.00861 }
\author{Michael K. Borregaard, Thomas J. Matthews & Karl I. Ugland}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create_octaves}}
}
\examples{
data(moths)
fit <- fitGambin(moths) 
plot(fit)
summary(fit)

# gambin parameters based on a standardized sample size of 1000 individuals
stand_fit <- replicate(20, fitGambin(moths, 1000)$Alpha) #may take a while on slower computers
print(c(mean = mean(stand_fit), sd = sd(stand_fit)))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


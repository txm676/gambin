% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_gambin.R
\name{summary.gambin}
\alias{summary.gambin}
\title{Summarising the results of a gambin model fit}
\usage{
\method{summary}{gambin}(object, confint = FALSE, n = 50, ...)
}
\arguments{
\item{object}{A gambin model fit object from \code{fit_abundances}}

\item{confint}{A logical argument specifying whether confidence intervals
should be calculated (via bootstrapping) for the parameters of gambin
models with more than 1 component (confidence intervals for 1 component
gambin models are calculated automatically)}

\item{n}{The number of bootstrap samples to use in generating the confidence
intervals (for multimodal gambin models)}

\item{\dots}{Further arguments to pass}
}
\value{
A list of class 'summary.gambin' with nine elements, containing
  useful information about the model fit.
}
\description{
S3 method for class 'gambin'. \code{summary.gambin} creates
  summary statistics for objects of class 'gambin'.The summary method
  generates more useful information (e.g. confidence intervals) for the user
  than the standard model fitting function. Another S3 method
  (\code{print.summary.gambin}; not documented) is used to print the output.
}
\details{
For the one-component gambin model the confidence interval for the
  alpha parameter is calculated automatically using an analytical solution.

  For gambin models with more than one component no analytical solution for
  deriving the confidence intervals is known. Instead, a bootstrapping
  procedure can be used (using the \code{confint} and \code{n} arguments) to
  generate confidence intervals around the alpha and max octave parameters.
  However, the process can be time-consuming, particularly for gambin models
  with more than two components. Thus, the default is that confidence
  intervals are not automatically calculated for gambin models with more than
  one component (i.e. \code{confint} == FALSE).

  In addition, it should be noted that in certain case the confidence
  intervals around the alpha parameters in multi-component gambin models can
  be quite wide. This is due to changes in the max octaves of the component
  distributions in the bootstrapped samples. It can be useful to make a plot
  (e.g. a dependency boxplot) of the n alpha values against the max octave
  values.
}
\examples{
 \dontrun{
data(moths)
fit = fit_abundances(moths)
summary(fit)
# multimodal gambin models with confidence intervals
biMod <- fit_abundances(moths, no_of_components = 2)
summary(biMod, confint = TRUE, n = 5) #large n takes a long time to run
}
}

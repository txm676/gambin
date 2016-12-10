#' @title Likelihood statistics for the GamBin model
#' @description Uses likelihood and information theoretical approaches to reveal 
#' the degree of fit of the GamBin model to empirical species abundance distributions.
#' @param object An object of type \code{gambin}
#' @param \dots Further arguments to pass to the function 
#' @return logLik returns an R object of type \code{logLik}. The other function return the numerical value of the statistic
#' @references Akaike, Hirotugu. "A new look at the statistical model identification." Automatic Control, 
#' IEEE Transactions on 19.6 (1974): 716-723.
#' @examples 
#' data(moths)
#' fit <- fitGambin(moths)
#' AIC(fit)
#' @importFrom stats AIC chisq.test coef confint logLik
#' @importFrom stats nobs optimise pgamma predict qchisq qgamma
#' @export
logLik.gambin <-
function(object, ...)
{
  return(object$logLik)
}

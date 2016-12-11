logLik_gamBin = function(alpha, mydata, maxoct) #the maxoct argument can be removed if we do not estimate maxoct
{
  if(missing(maxoct)) maxoct <- max(mydata$octave)
  dgamb <- dgambin(alpha, maxoct)
  exponent <- mydata$species # this line and the next can be removed if we do not estimate maxoctave
  if(length(exponent) < length(dgamb)) exponent[(length(exponent)+1):length(dgamb)] <- 0
  lik_dist <- dgamb^exponent
  logLik <- sum(-log(lik_dist))
  if(logLik == Inf)
    return(999999999999)
  logLik
}

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
logLik.gambin = function(object, ...)
{
  return(object$logLik)
}

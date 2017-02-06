


ll_w = function(alpha, w, maxoctave, values, freq) {
  (res = vapply(seq_along(alpha), 
                function(i) w[i]*dgambin(values, alpha[i], maxoctave, log=FALSE), 
                FUN.VALUE = numeric(length(values))))
  log_prob = log(rowSums(res))   
  sum(freq*log_prob)
}

## This function extracts the parameters from par
## Then sanity checks before passing on to ll_w
ll_optim = function(par, maxoctave, values, freq) {
  len_p = (length(par) + 1)/2
  alpha=par[1:len_p]; w=par[(len_p+1):(2*len_p)]
  w[len_p] = 1 - sum(w, na.rm=TRUE)
  if(any(alpha < 0) || any(w < 0) || any(w > 1)) return(Inf)
  #w = w/sum(w)
  -ll_w(alpha, w, maxoctave, values, freq)  
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
#' fit = fit_abundances(moths)
#' AIC(fit)
#' @importFrom stats AIC chisq.test coef confint logLik
#' @importFrom stats nobs optimise pgamma predict qchisq qgamma
#' @export
logLik.gambin = function(object, ...)
{
  return(object$logLik)
}

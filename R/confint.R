logLik_gamBin = function(alpha, x) {
  dgamb = dgambin(x$octave, alpha = alpha, maxoctave = max(x$octave), 
                  w = 1, log = TRUE)
  sum(x$species * dgamb)
}

est_confint  = function(est_alpha, est_likelihood, mydata, level) {
  
  conf_logLik = function(alpha, mydata, est_likelihood, level)
  {
    # the likelihood ratio test for confidence intervals (from "Beyond Traditional Statistical Measures")
    # -2 LL(pi)/ll(hat{pi}) < chi{1}
    return(abs(-logLik_gamBin(alpha, mydata) + est_likelihood - qchisq(level,1)/2))
  }

  lower = optimise(conf_logLik, interval = c(0, est_alpha),
                   mydata = mydata, est_likelihood = est_likelihood, level = level)$minimum
  
  upper_alpha = 4*((est_alpha - lower) + est_alpha) # Crude upper estimate
  higher = optimise(conf_logLik, interval = c(est_alpha, upper_alpha),
                    mydata = mydata, est_likelihood = est_likelihood, level = level)$minimum
  
  return(c(lower, higher))
}

#' @export 
confint.gambin = function(object, parm = "alpha", level = 0.95, ...)
{
  if(!tolower(parm) == "alpha") stop("Only the alpha parameter has confidence intervals", call. = FALSE)
  est_confint(est_alpha = object$alpha,
              est_likelihood = object$logLik,
              mydata = object$Data, 
              level = level)
}

#' @export
summary.gambin = function(object, ...)
{
  res = list()
  res$Data = object$Data
  res$Dataname = object$Dataname
  
  res$alpha = object$alpha
  res$octaves = object$octaves
  res$coefficients = object$coefficients
  res$residuals = object$Data$species - object$fitted.values
  
  res$ConfInt95 = if(length(res$alpha) > 1) NA else confint(object)
  res$logLik = object$logLik
  
  chiprobs = object$fitted.values/sum(object$fitted.values)
  suppressWarnings(res$ChiSq <- chisq.test(object$Data$species, p = chiprobs))
  
  attr(res, "nobs") = nrow(res$Data)
  class(res) = "summary.gambin"
  return(res)
}

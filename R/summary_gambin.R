#' @title Summarising the results of a gambin model fit
#' @description S3 method for class 'gambin'. \code{summary.gambin} creates
#'   summary statistics for objects of class 'gambin'.The summary method
#'   generates more useful information (e.g. confidence intervals) for the user
#'   than the standard model fitting function. Another S3 method
#'   (\code{print.summary.gambin}; not documented) is used to print the output.
#' @param object A gambin model fit object from \code{fit_abundances}
#' @param confint A logical argument specifying whether confidence intervals
#'   should be calculated (via bootstrapping) for the parameters of gambin
#'   models with more than 1 component (confidence intervals for 1 component
#'   gambin models are calculated automatically)
#' @param n The number of bootstrap samples to use in generating the confidence
#'   intervals (for multimodal gambin models)
#' @param \dots Further arguments to pass
#' @details For the one-component gambin model the confidence interval for the
#'   alpha parameter is calculated automatically using an analytical solution.
#'
#'   For gambin models with more than one component no analytical solution for
#'   deriving the confidence intervals is known. Instead, a bootstrapping
#'   procedure can be used (using the \code{confint} and \code{n} arguments) to
#'   generate confidence intervals around the alpha and max octave parameters.
#'   However, the process can be time-consuming, particularly for gambin models
#'   with more than two components. Thus, the default is that confidence
#'   intervals are not automatically calculated for gambin models with more than
#'   one component (i.e. \code{confint} == FALSE).
#'
#'   In addition, it should be noted that in certain case the confidence
#'   intervals around the alpha parameters in multi-component gambin models can
#'   be quite wide. This is due to changes in the max octaves of the component
#'   distributions in the bootstrapped samples. It can be useful to make a plot
#'   (e.g. a dependency boxplot) of the n alpha values against the max octave
#'   values.
#' @return A list of class 'summary.gambin' with nine elements, containing
#'   useful information about the model fit.
#' @importFrom stats chisq.test
#' @examples
#' data(moths)
#' fit = fit_abundances(moths)
#' summary(fit)
#' # multimodal gambin models with confidence intervals
#' biMod <- fit_abundances(moths, no_of_components = 2)
#' summary(biMod, confint = TRUE, n = 5) #large n takes a long time to run
#' @export
#' 
summary.gambin = function(object, confint = FALSE, n = 50, ...)
{
  res = list()
  res$Data = object$Data
  res$Dataname = object$Dataname
  
  res$alpha = object$alpha
  res$octaves = object$octaves
  res$coefficients = object$coefficients
  res$residuals = object$Data$species - object$fitted.values
  
  res$ConfInt95 = if(length(res$alpha) == 1){
    confint(object)
  } else {
    if (!confint){
      NA
    } else {
      suppressMessages(confint(object, n = n, no_of_components = length(res$alpha), cores = object$cores))
    }
  }
  res$logLik = object$logLik
  
  chiprobs = object$fitted.values/sum(object$fitted.values)
  suppressWarnings(res$ChiSq <- chisq.test(object$Data$species, p = chiprobs))
  
  attr(res, "nobs") = nrow(res$Data)
  class(res) = "summary.gambin"
  return(res)
}

logLik_gamBin = function(alpha, x) {
  dgamb = dgambin(x$octave, alpha = alpha, maxoctave = max(x$octave), 
                  w = 1, log = TRUE)
  sum(x$species * dgamb)
}

est_confint  = function(est_alpha, est_likelihood, sample_data, level) {
  
  conf_logLik = function(alpha, sample_data, est_likelihood, level)
  {
    # the likelihood ratio test for confidence intervals (from "Beyond Traditional Statistical Measures")
    # -2 LL(pi)/ll(hat{pi}) < chi{1}
    return(abs(-logLik_gamBin(alpha, sample_data) + est_likelihood - qchisq(level,1)/2))
  }

  lower = optimise(conf_logLik, interval = c(0, est_alpha),
                   sample_data = sample_data, est_likelihood = est_likelihood, level = level)$minimum
  
  upper_alpha = 4*((est_alpha - lower) + est_alpha) # Crude upper estimate
  higher = optimise(conf_logLik, interval = c(est_alpha, upper_alpha),
                    sample_data = sample_data, est_likelihood = est_likelihood, 
                    level = level)$minimum
  
  return(c(lower, higher))
}

# parm for S3 generic consistency
#' @export 
confint.gambin = function(object, parm = "alpha", level = 0.95, n = n, 
                          no_of_components = no_of_components, cores = cores, ...)
{
  if (length(object$alpha) == 1){
    if(!tolower(parm) == "alpha") stop("Only the alpha parameter has confidence intervals", call. = FALSE)
    est_confint(est_alpha = object$alpha,
                est_likelihood = object$logLik,
                sample_data = object$Data, 
                level = level)
  } else {
    ab = object$Data
    ci = replicate(n, bs_confint(ab, no_of_components = no_of_components , cores = cores))
    d = round(apply(ci, 1, function(x) quantile(x, c(0.025, 0.975))), 3)
    d = as.data.frame(d)
    n = ncol(d)
    s1 = seq(1, (n/2), 1)
    colnames(d)[1:(n/2)] = paste("alpha", s1, sep = "")
    colnames(d)[((n/2) + 1):n] = paste("max_octave", s1, sep = "")
    d
  }
}



##bootstrappoed CIs for models with > 1 component
bs_confint = function(abundances, no_of_components = 1, cores = 1) {
  x = sample(abundances$octave, 
             size = sum(abundances$species), 
             prob = abundances$species, replace = TRUE)
  x = table(x)
  freq = as.vector(x)
  values = as.numeric(as.character(names(x)))
  abundances = data.frame(octave=values, species = freq)
  
  alpha.boot <- fit_abundances(abundances, no_of_components = no_of_components, 
                               cores = cores)
  return(c(alpha.boot$alpha, alpha.boot$octaves))
}


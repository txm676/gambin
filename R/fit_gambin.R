core_message = function(cores) {
  if(cores == 1L) {
    message("Using 1 core. Your machine has ", parallel::detectCores(), " available.")
  } else {
    message("Using ", cores, " cores. Your machine has ", parallel::detectCores(), " cores available.")
  }
}


#' @title Fit a unimodal or multimodal gambin model to a species abundance distribution
#' @description Uses maximum likelihood methods to fit the GamBin model (with a given number of modes) to binned
#' species abundances. To control for the effect of sample size, the abundances
#' may be subsampled prior to fitting.
#' @param abundances Either a vector of abundances of all species in the sample/community; or the result of \code{create_octaves}
#' @param subsample The number of individuals to sample from the community before fitting the GamBin model.
#' If subsample == 0 the entire community is used
#' @param no_of_components Number of components (i.e. modes) to fit.The default (no_of_components == 1) fits the standard
#' unimodal gambin model.
#' @param \dots Further arguments to pass to \code{barplot}
#' @param cores No of cores to use when fitting. Use \code{parallel::detectCores()} to
#' detect the number of cores on your machine. 
#' @details The gambin distribution is fit to the number of species in abundance octaves,
#' as specified by the \code{create_octaves} function. Because the shape of species abundance
#' distributions depend on sample size, abundances of different communities should be compared
#' on equally large samples. The sample size can be set by the \code{subsample} parameter.
#' To estimate \code{alpha} from a standardised sample, the function must be run several
#' times; see the examples. The \code{no_of_components} parameter enables mutlimodal gambin
#' distributions to be fitted. For example, setting \code{no_of_components} equal to 2, the bimodal
#' gambin model is fitted. When a multimodal gambin model is fitted (with g modes), the return values are the alpha
#' parameters of the g different component distributions, the max octave values for the g component distributions 
#' (as the max octave values for the g-1 component distributions are allowed to vary), and the and the weight parameter(s) 
#' which denote the fraction of objects within each g component distribution. When fitting multimodal gambin models
#' (particuarly on large datasets), the optimisation algorithm can be slow. In such cases, the process
#' can be speeded up by using the \code{cores} parameter to enable parallel computing.
#' 
#' The \code{plot} method creates a barplot showing the observed
#' number of species in octaves, with the fitted GamBin distribution shown as black dots.
#' @return The \code{fit_abunbances} function returns an object of class \code{gambin}, with the \code{alpha},
#' \code{w}
#' and \code{MaxOctave} parameters of the GamBin mixture distribution,
#' the likelihood of the fit, and the empirical distribution over octaves.
#' @importFrom stats optim
#' @examples
#' data(moths)
#' fit = fit_abundances(moths)
#' barplot(fit)
#' lines(fit, col=2)
#' summary(fit)
#' # gambin parameters based on a standardized sample size of 1000 individuals
#' stand_fit <- replicate(20, fit_abundances(moths, 1000)$alpha) #may take a while on slower computers
#' print(c(mean = mean(stand_fit), sd = sd(stand_fit)))
#' # a bimodal gambin model
#' biMod <- fit_abundances(moths, no_of_components = 2)
#' @export
fit_abundances <- function(abundances, subsample = 0, no_of_components = 1, cores = 1)
{

  # test for NA's in the data (and that it is a numeric vector)
  Dataname <- deparse(substitute(abundances))

  if(is.vector(abundances) && is.numeric(abundances)) {
    mydata = create_octaves(abundances, subsample)
  } else  if(is.data.frame(abundances)) {
    names(abundances) <- tolower(names(abundances))
    if(!("species" %in% names(abundances) && "octave" %in% names(abundances))) {
      stop("abundances must be a numeric vector or a data.frame created by create_octaves")
    }
    mydata = abundances[c("octave", "species")]
  } else {
    stop("abundances must be a numeric vector or a data.frame created by create_octaves")
  }

  if(no_of_components == 1) {
    val = optim(par = 1, fn=ll_optim, maxoctave = max(mydata$octave),
                values = mydata$octave, freq = mydata$species, method = "Brent",
                lower = 0, upper = 100)
    val$octaves = max(mydata$octave)
    ## Optim minimises so convert to negative
    logLik = -val$value
  } else {
    core_message(cores)
    alpha = rep(1, no_of_components)
    w = rep(1/length(alpha), length(alpha) - 1)
    val = estimate_parameters(par = c(alpha, w), values = mydata$octave, 
                              freq = mydata$species, cores = cores)
    logLik = -val$value
  }
  res = list()
  res$alpha = val$par[1:no_of_components]
  res$w = val$par[(no_of_components+1):(2*no_of_components)]
  res$w[length(res$w)] = 1-sum(res$w, na.rm = TRUE)
  res$octaves = val$octaves
  res$convergence = val$convergence
  
  ## -1: weights have to sum to 1
  ## -1: one of the octaves is equal to max(data)
  attr(logLik, "df") = no_of_components*2 + (no_of_components - 1) 
  attr(logLik, "nobs") = nrow(mydata)
  class(logLik) = "logLik"
  res$logLik = logLik

  res$fitted.values =  dgambin(mydata$octave,
                               alpha = res$alpha, w=res$w,
                               maxoctave = res$octaves) * sum(mydata$species)
  res$Data = mydata
  if (subsample == 0)
    res$Dataname = Dataname else
      res$Dataname = paste(subsample, "individuals sampled from", Dataname)
  if(length(res$Dataname) > 1) res$Dataname = "unnamed abundance vector"

  res$max_octaves = res$max_octaves
  res$coefficients = c(alpha = res$alpha, w = res$w, max_octave = res$octaves)
  
  #get fitted values of each component distribution and the peak value of each for deconstruction function
  if (no_of_components > 1) {
    fvals <- mapply(function(x, y) dgambin_single(res$Data$octave, x, y) * 
                                    sum(mydata$species), x = res$alpha, y = res$octaves)
    res$peak_octave_individual <- apply(fvals, 2, which.max) - 1 #minus one as octaves start at 0
  }
  
  attr(res, "nobs") = nrow(mydata)
  class(res) = "gambin"
  res
}

#' @rdname fit_abundances
#' @export
fitGambin = function(abundances, subsample = 0) {
  .Deprecated("fit_abundances")
  fit_abundances(abundances = abundances, subsample)
}












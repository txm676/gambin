core_message = function(cores) {
  if(cores == 1L) {
    message("Using 1 core. Your machine has ", parallel::detectCores(), " available.")
  } else {
    message("Using ", cores, " cores. Your machine has ", parallel::detectCores(), " cores available.")
  }
}


#' @title Fit the GamBin model to a species abundance distribution
#' @description Uses maximum likelihood methods to fit the GamBin model to binned
#' species abundances. To control for the effect of sample size, the abundances
#' may be subsampled prior to fitting.
#' @param abundances Either a vector of abundances of all species in the sample/community; or the result of \code{create_octaves}
#' @param subsample The number of individuals to sample from the community before fitting the GamBin model.
#' If subsample == 0 the entire community is used
#' @param no_of_components Number of components to fit.
#' @param \dots Further arguments to pass to \code{barplot}
#' @param cores No of cores to use when fitting. Use \code{parallel::detectCores()} to
#' detect the number of cores on your machine. 
#' @details The gambin distribution is fit to the number of species in abundance octaves,
#' as specified by the \code{create_octaves} function. Because the shape of species abundance
#' distributions depend on sample size, abundances of different communities should be compared
#' on equally large samples. The sample size can be set by the \code{subsample} parameter.
#' To estimate \code{alpha} from a standardised sample, the function must be run several
#' times; see the examples. The \code{plot} method creates a barplot showing the observed
#' number of species in octaves, with the fitted GamBin distribution shown as black dots.
#' @return The \code{fit_abunbances} function returns an object of class \code{gambin}, with the \code{alpha},
#' \code{w}
#' and \code{MaxOctave} parameters of the GamBin mixture distribution,
#' the likelihood of the fit, and the empirical distribution over octaves.
#' @importFrom stats optim
#' @examples
#' data(moths)
#' fit = fit_abundances(create_octaves(moths))
#' barplot(fit)
#' lines(fit, col=2)
#' summary(fit)
#' # gambin parameters based on a standardized sample size of 1000 individuals
#' stand_fit <- replicate(20, fit_abundances(moths, 1000)$alpha) #may take a while on slower computers
#' print(c(mean = mean(stand_fit), sd = sd(stand_fit)))
#' @export
fit_abundances <- function(abundances, subsample = 0, no_of_components = 1, cores = 1)
{

  # test for NA's in the data (and that it is a numeric vector)
  Dataname <- deparse(substitute(abundances))

  if(is.vector(abundances) && is.numeric(abundances)) {
    warning("Calling create_octaves on input. This feature is Deprecated. Please call
            create_octaves explicitly.")
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
    val = optim(par = 1, fn=ll_optim, maxoct = max(mydata$octave),
                values = mydata$octave, freq = mydata$species, method = "Brent",
                lower = 0, upper = 100)
    val$octaves = max(mydata$octave)
  } else {
    core_message(cores)
    alpha = rep(1, no_of_components)
    w = rep(1/length(alpha), length(alpha) - 1)
    val = estimate_parameters(par = c(alpha, w), values = mydata$octave, 
                              freq = mydata$species, cores = cores)
  }
  res = list()
  res$alpha = val$par[1:no_of_components]
  res$w = val$par[(no_of_components+1):(2*no_of_components)]
  res$w[length(res$w)] = 1-sum(res$w, na.rm = TRUE)
  res$octaves = val$octaves
  res$convergence = val$convergence
  
  logLik = -val$value
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












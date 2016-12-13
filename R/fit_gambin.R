#' @title Fit the GamBin model to a species abundance distribution
#' @description Uses maximum likelihood methods to fit the GamBin model to binned 
#' species abundances. To control for the effect of sample size, the abundances 
#' may be subsampled prior to fitting. 
#' @param abundances Either a vector of abundances of all species in the sample/community; or the result of \code{create_octaves} 
#' @param subsample The number of individuals to sample from the community before fitting the GamBin model. If subsample == 0 the entire community is used
#' @param object a \code{gambin} object created by \code{fitGambin}
#' @param x A \code{gambin} object created by \code{fitGambin}
#' @param parm The parameter to calculate confidence intervals from. Only alpha is implemented
#' @param level The significance level of the confidence intervals 
#' @param barcol The colour of the bars illustrating the empirical abundance of species in octaves 
#' @param barwidth The width of the bars illustrating the empirical abundance of species in octaves 
#' @param cex.dots The size of the dots illustrating the fitted abundance of species in octaves 
#' @param dotpch The point character of the dots illustrating the fitted abundance of species in octaves 
#' @param dotcol The colour of the dots illustrating the fitted abundance of species in octaves 
#' @param line Should the dots be connected with a line? 
#' @param lwd The width of the line connecting dots 
#' @param linecol The colour of the line connecting dots 
#' @param \dots Further arguments to pass to \code{barplot}
#' @details The gambin distribution is fit to the number of species in abundance octaves, 
#' as specified by the \code{create_octaves} function. Because the shape of species abundance 
#' distributions depend on sample size, abundances of different communities should be compared 
#' on equally large samples. The sample size can be set by the \code{subsample} parameter. 
#' To estimate \code{alpha} from a standardised sample, the function must be run several 
#' times; see the examples. The \code{plot} method creates a barplot showing the observed 
#' number of species in octaves, with the fitted GamBin distribution shown as black dots.
#' @return The \code{fitGambin} function returns an object of class \code{gambin}, with the \code{alpha} 
#' and \code{MaxOctave} parameters of the GamBin distribution, the likelihood of the fit, and the empirical distribution over octaves.
#' @examples 
#' data(moths)
#' fit <- fitGambin(moths) 
#' plot(fit)
#' summary(fit)
#' # gambin parameters based on a standardized sample size of 1000 individuals
#' stand_fit <- replicate(20, fitGambin(moths, 1000)$Alpha) #may take a while on slower computers
#' print(c(mean = mean(stand_fit), sd = sd(stand_fit)))
#' @export
fitGambin <- function(abundances, subsample = 0)
{
  # for the GamBin function, all the abundances are binned into octaves (on log base 2)
  
  # test for NA's in the data (and that it is a numeric vector)
  Dataname <- deparse(substitute(abundances))
  
  if(is.vector(abundances) && is.numeric(abundances)) {
    mydata <- create_octaves(abundances, subsample)
  } else  if(is.data.frame(abundances)) {
    names(abundances) <- tolower(names(abundances))
    if(!("species" %in% names(abundances) && "octave" %in% names(abundances))) {
      stop("abundances must be a numeric vector or a data.frame created by create_octaves")	
    }
    mydata <- abundances[c("octave", "species")]
  } else {
    stop("abundances must be a numeric vector or a data.frame created by create_octaves")	
  }
  
  val <- optimise(logLik_gamBin, interval = c(0, 15), mydata = mydata)
  res <- list()
  
  res$Alpha <- val[[1]]
  if(res$Alpha == 15) 
  {
    warning("Alpha could not be estimated")
    res$Alpha <- Inf
  }
  
  logLik <- -val[[2]]
  attr(logLik, "df") <- 1
  attr(logLik, "nobs") <- nrow(mydata)
  class(logLik) <- "logLik"
  res$logLik <- logLik
  
  res$fitted.values <- gambin_exp(res$Alpha, max(mydata$octave), sum(mydata$species))
  res$Data <- mydata
  if (subsample == 0)
    res$Dataname <- Dataname else
      res$Dataname <- paste(subsample, "individuals sampled from", Dataname)
  if(length(res$Dataname) > 1) res$Dataname <- "unnamed abundance vector"
  
  res$MaxOctave <- max(mydata$octave)
  
  res$coefficients <- c(Alpha = res$Alpha, MaxOctave = res$MaxOctave)
  
  attr(res, "nobs") <- nrow(mydata)
  class(res) <- "gambin"
  res
}

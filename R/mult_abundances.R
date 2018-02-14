#' @title Fit a unimodal gambin model to multiple species abundance distributions
#' @description Fits the unimodal gambin model to the SADs from multiple sites and
#' returns the standardised and unstandardised alpha values. 
#' @param mult Either a matrix, dataframe or list containing the species abundance data of a set of sites.
#' In the case of a matrix or dataframe, a given column contains the abundance data for a given site (i.e. columns 
#' are sites and rows are species; each cell is the abundance of a given species in a given site). In the case
#' of a list, each element in the list contains the abundance data (i.e. a vector of abundances) for a given site.
#' @param N The number of times to subsample the abundance data in each site to calculate mean standardised alpha.
#' @param subsample The number of individuals to sample from each site before fitting the gambin model. The default is
#' subsample = NULL, in which case subsample is set to equal the number of individuals in the site with the 
#' fewest individuals.
#' @param r The number of decimal points to round the returned alpha values to (default is r = 3)
#' @details Because the alpha parameter of the gambin model is dependent on sample size, when comparing the alpha 
#' values between sites it can be useful to first standardise the number of individuals in all sites. By default, 
#' the \code{mult_abundances} function calculates the total number of individuals in each site and selects the 
#' minimum value for standardising. This minimum number of individuals is then sampled from each site and the 
#' gambin model fitted to this subsample using \code{fit_abundances}, and the alpha value stored. This process is 
#' then repeated \code{N} times and the mean alpha value is calculated for each site. The number of individuals 
#' to be subsampled can be manually set using the \code{subsample} argument. The function returns a list, in which 
#' the first two elements are the mean standardised alpha values for each site, and the raw unstandardized alpha 
#' values for each site, respectively. The full set of \code{N} alpha values and X2 P-values for each site are 
#' also returned. 

#' As an input, the SAD data can be in the form of a matrix or dataframe, or a list. A matrix/dataframe is only for 
#' when each site (column) has abundance data for the same set of species (rows). For example, an abundance matrix of 
#' bird species on a set of islands in an archipelago. A list should be used in cases where the number of species 
#' varies between sites; for example, when comparing the SADs of samples from different countries. In this case, 
#' each element of the list contains an abundance vector from a given site. 
#' 
#' At present, the \code{mult_abundances} function only fits the unimodal gambin model. 
#' @examples
#' #simulate a matrix containing the SAD data for 20 sites (50 sp. in each)
#' mult <- matrix(0, nrow = 50, ncol = 20)
#' mult <- apply(mult, 2, function(x) ceiling(rlnorm(length(x), 0, 2.5)))
#' 
#' #run the mult_abundances function and view the alpha values
#' mm <- mult_abundances(mult, N = 20, subsample = NULL)
#' mm[1:2]
#' plot(mm$Mean.Stan.Alpha, mm$Unstan.Alpha)
#' 
#' #simulate a list containing the SAD of 20 sites (with varying numbers of sp.)
#' mult2 <- vector("list", length = 20)
#' for (i in 1:ncol(mult)){
#'  dum <- sample(mult[, i], replace = TRUE)
#'  rm <- round(runif(1, 0, 5), 0)
#'  if (rm > 0){
#'    rm2 <- sample(1:length(dum), rm, replace = FALSE)
#'    dum <- dum[-rm2]
#'  }
#'  mult2[[i]] <- dum
#' }
#' 
#' #run the mult_abundances function on the list
#' mm2 <- mult_abundances(mult2, N = 20, subsample = NULL)
#' mm2[1:2] 
#' @export


mult_abundances <- function(mult, N = 100, subsample = NULL,  r = 3){
  
  if (!(is.list(mult) || is.matrix(mult) || is.data.frame(mult))) stop("Data are not in correct format")
  if (is.data.frame(mult)) mult <- as.matrix(mult)
  if (is.list(mult)) if (any(vapply(mult, anyNA, logical(1)))) stop("NAs present in data")
  if (is.matrix(mult)) if (anyNA(mult)) stop("NAs present in data")
  
  #if mult was either a matrix or df (df changed to matrix above)
  if (is.matrix(mult)){
    if (is.null(subsample)){
      minAbun <- min(colSums(mult))[1]
      message ("\n", "Minimum number of individuals in a site = ", minAbun,"\n", 
               "Subsampling ", minAbun, " individuals from each site", "\n", sep = " ")
      if (minAbun < 50) warning("Minimum abundance is below 50: this is very small!", "\n")
      multfit <- replicate(N, apply(mult, 2, fit_abundances, subsample = minAbun, 
                                    no_of_components = 1, cores = 1))
    }
    
    if (!is.null(subsample)){
      minAbun2 <- min(colSums(mult))[1]
      if (subsample > minAbun2) stop("User supplied subsample value exceeds at least one site's total abundance")
      minAbun <- subsample
      message ("\n", "Minimum number of individuals in a site = ", minAbun2,"\n", 
               "Subsampling ", subsample, " individuals from each site", "\n", sep = " ")
      if (minAbun < 50) warning("Minimum abundance is below 50: this is very small!", "\n")
      multfit <- replicate(N, apply(mult, 2, fit_abundances, subsample = minAbun, 
                                    no_of_components = 1, cores = 1))
    }
    
    #standardised alpha
    alpV <- round(vapply(multfit, function(x) x$alpha, FUN.VALUE = numeric(1)), r)
    alpMat <- matrix(alpV, nrow = N, ncol = ncol(mult), byrow = TRUE)
    meanAlpha <- round(colMeans(alpMat), r)
    meanAlpha <- as.vector(meanAlpha)
    
    #unstandaridised alpha
    unAlp <- round(apply(mult, 2, function(x) (fit_abundances(x))$alpha), r)
    
    
    #X2 results
    x2v <- round(vapply(multfit, function(x) summary(x)$ChiSq$p.value, FUN.VALUE = numeric(1)), r)
    x2Mat <- matrix(x2v, nrow = N, ncol = ncol(mult), byrow = TRUE)
    if (any(x2Mat < 0.05)) warning("Some X2 P-value less than 0.05")
    
    
    resList <- list("Mean.Stan.Alpha" = meanAlpha, "Unstan.Alpha" = unAlp, "All.Alpha" = alpMat, "X2" = x2Mat)
  }
  
  #if mult is a list
  if (is.list(mult)){
    
    if (is.null(subsample)){
      minAbun <- min(vapply(mult, sum, FUN.VALUE = numeric(1)))[1]
      message ("\n", "Minimum number of individuals in a site = ", minAbun,"\n", 
               "Subsampling ", minAbun, " individuals from each site", "\n", sep = " ")
      if (minAbun < 50) warning("Minimum abundance is below 50: this is very small!", "\n")
      multfit <- replicate(N, lapply(mult, fit_abundances, subsample = minAbun, 
                                     no_of_components = 1, cores = 1))
    }
    
    if (!is.null(subsample)){
      minAbun2 <- min(vapply(mult, sum, FUN.VALUE = numeric(1)))[1]
      if (subsample > minAbun2) stop("User supplied subsample value exceeds at least one site's total abundance")
      minAbun <- subsample
      message ("\n", "Minimum number of individuals in a site = ", minAbun2,"\n", 
               "Subsampling ", subsample, " individuals from each site", "\n", sep = " ")
      if (minAbun < 50) warning("Minimum abundance is below 50: this is very small!", "\n")
      multfit <- replicate(N, lapply(mult, fit_abundances, subsample = minAbun, 
                                     no_of_components = 1, cores = 1))
    }
    
    
    #standardised alpha
    alpV <- round(vapply(multfit, function(x) x$alpha, FUN.VALUE = numeric(1)), r)
    alpMat <- matrix(alpV, nrow = N, ncol = length(mult), byrow = TRUE)
    meanAlpha <- round(colMeans(alpMat), r)
    meanAlpha <- as.vector(meanAlpha)
    
    #unstandaridised alpha
    unAlp <- round(vapply(mult, function(x) (fit_abundances(x))$alpha, FUN.VALUE = numeric(1)), r)
    
    
    #X2 results
    x2v <- round(vapply(multfit, function(x) summary(x)$ChiSq$p.value, FUN.VALUE = numeric(1)), r)
    x2Mat <- matrix(x2v, nrow = N, ncol = length(mult), byrow = TRUE)
    if (any(x2Mat < 0.05)) warning("Some X2 P-value less than 0.05")
    
    resList <- list("Mean.Stan.Alpha" = meanAlpha, "Unstan.Alpha" = unAlp, "All.Alpha" = alpMat, "X2" = x2Mat)
    
  }
  return(resList)
}

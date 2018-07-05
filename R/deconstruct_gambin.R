

#internal function similar to 'table' but which includes the frequencies of categories
#that are not in x but are in nams. This is needed for cases when a octave (x) is
#missing a category that is found in the dataset but just not that octave.
#nams is a complete vector of unique categories in the dataset
table_int <- function(x, nams){
  x1 <- table(x)
  if (length(x1) == length(nams)) {
    x1 <- as.data.frame(x1)
    return(x1)}
  n1 <- names(x1)
  n2 <- setdiff(nams, n1)
  x1 <- as.data.frame(x1)
  n3 <- data.frame(n2, 0)
  colnames(n3) <- colnames(x1)
  x2 <- rbind(x1, n3)
  return(x2)
}



#print function

print.deconstruct <- function(x){
  cat("\n")
  cat("Multimodal gambin modal octaves deconstruction analysis\n\n")
  cat(paste("Number of modes in gambin fit: ", length(x$Peak_locations), "\n", sep = ""))
  cat(paste("Peak octave locations: ", paste(x$Peak_locations,  collapse = " "),"\n\n", sep=""))
  if (length(x) == 2) cat("No species classification data provided\n")
  if (length(x) == 3) {
    cat("Species classification summary table:\n\n")
    print(x$Summary_table)
  }
}



#' @title Deconstruct a multimodal gambin model fit
#' @description Deconstruct a multimodal gambin model fit by locating the modal
#'   octaves and (if species classification data are provided) determining the
#'   proportion of different types of species in each octave.
#' @param fit A gambin model fit where the number of components is greater than
#'   one (see \code{\link{fit_abundances}}).
#' @param dat A matrix or dataframe with at least two columns, including the
#'   abundance data used to fit the multimodal gambin model and the species
#'   names. An optional third column can be provided that contains species
#'   classification data.
#' @param peak_val A vector of of modal octave values. If \code{peak_val =
#'   NULL}, the modal octave values are taken from the model fit object.
#' @param abundances The name of the column in \code{dat} that contains the
#'   abundance data (default = "abundance").
#' @param species The name of the column in \code{dat} that contains the species
#'   names (default = "species").
#' @param categ Either NULL if no species classification data are provided, or
#'   the name of the column in \code{dat} that contains the species
#'   classification data.
#' @param plot_modes A logical argument specifying whether a barplot of the
#'   model fit with highlighted octaves should be generated. If \code{categ =
#'   FALSE} a barplot is produced wherby just the modal octaves are highlighted
#'   in red. If \code{categ = TRUE} a barplot is produced whereby the bar for
#'   each octave is split into n parts, where n equals the number of species
#'   categories.
#' @param col.statu A vector of colours (of length n) for the split barplot,
#'   where n equals the number of species categories.
#' @details The function enables greater exploration of a multimodal gambin
#'   model fit. If no species classification data are available (i.e.
#'   \code{categ = NULL}) the function returns the modal octaves of the
#'   n-component distributions and the names of the species located in each
#'   octave. If \code{plot_modes = TRUE} a plot is returned with the modal
#'   octaves highlighted in red. If species classification data are provided the
#'   function also returns a summary table with the number of each species
#'   category in each octave provided. The user can then use these data to run
#'   different tests to test whether, for example, the number of species in each
#'   category in the modal octaves is signficantly different than expected by
#'   chance. If \code{plot_modes = TRUE} a split barplot is returned whereby
#'   each bar (representing an octave) is split into the n species categories.
#'
#'   Species classification data should be of type character (e.g. native or
#'   invasive).
#'
#'   Occasionally, some of the component distributions in a multimodal gambin
#'   model fit have the same modal octave; this is more common when fitting the
#'   3-component model. When this occurs a warning is produced, but it is not a
#'   substantive issue.
#' @return An object of class \code{deconstruct}. The object
#'   is a list with either two or three elements. If \code{categ = NULL}, the
#'   list has two elements: 1) 'Peak_locations', which contains the modal octave
#'   values, and 2) 'Species_per_octave', which is a list where each element
#'   contains the species names in an octave. If \code{categ != NULL}, the
#'   returned object has a third element: 3) 'Summary_table', which contains a
#'   dataframe (frequency table) with the numbers of each category of species in
#'   each octave.
#' @author Thomas J. Matthews & Francois Rigal
#' @importFrom graphics barplot legend points 
#' @importFrom grDevices rainbow
#' @examples
#' data(categ)
#' fits2 = fit_abundances(categ$abundances, no_of_components = 2)
#' #without species classification data
#' deconstruct_modes(fits2, dat = categ, peak_val = NULL, abundances = "abundances", 
#' species = "species", categ = NULL, plot_modes = TRUE)
#' #with species classification data
#' deconstruct_modes(fits2, dat = categ, categ = "status", col.statu = c("green", "red", "blue"))
#' #manually choose modal octaves
#' deconstruct_modes(fits2, dat = categ, peak_val = c(0,1))
#' @export

deconstruct_modes <- function(fit, dat, peak_val = NULL, abundances = "abundances", species = "species", categ = NULL,
                              plot_modes = TRUE, col.statu=NULL){
  
  if (class(fit) != "gambin") stop("fit is not a gambin object")
  if (!(is.data.frame(dat) || is.matrix(dat))) stop ("dat should be a dataframe or matrix")
  if (is.matrix(dat)) dat <- as.data.frame(dat)
  if (length(unique(fit$peak_octave_individual)) < length(fit$peak_octave_individual)) {
    warning("Two or more component distributions in the model fit have the same modal octave (see documentation)")
  }
  #if user does not provide peaks, take from the model fit (new code implemented in fit_abundanceS)  
  if (is.null(peak_val)) {
    #need to use sort as occasionally the peaks are not increasing order
    peaks <- sort(fit$peak_octave_individual, decreasing = FALSE)
  } else{
    peaks <- sort(peak_val, decreasing = FALSE)
  }
  
  oct_lim_min <- floor(2^(fit$Data$octave))
  oct_lim_max <- (oct_lim_min * 2) - 1
  oct_lim <- data.frame("min" = oct_lim_min, "max" = oct_lim_max)
  
  #create a list where each element contains a filtered 'dat' containing the species from a modal octave
  dat_peaks <- apply(oct_lim, 1, function(x) dat[dat[[abundances]] >= x[1] & dat[[abundances]] <= x[2],])
  
  #list of just the species names in each modal octave
  HE <- lapply(dat_peaks, function(x) as.vector(x[[species]]))
  names(HE) <- paste("Octave", fit$Data$octave, sep = "_")
  
  if (!is.null(categ)){
    #list where each element is the categories inside each  octave
    var_peaks <- lapply(dat_peaks, function(x) x[[categ]])
    nams <- unique(unlist(dat[[categ]]))#get all possible categories across the whole dataset
    
    #return a list where each element is a table (data frame) with the frequency of each category
    #in each  octave (including categories with 0 species, i.e. in the whole dataset but not a given mode)
    vpt1 <- lapply(var_peaks, table_int, nams = nams)
    ##create dataframes with number and proportion of each category in each mode
    #number
    vptD <- lapply(vpt1, function(x) { #turn category names from factor to character
      x$x <- as.character(x$x)
      x})
    vpt1 <- lapply(vptD, function(x) x[order(x$x),])#sort each element df by category name (alphabetical)
    #convert vpt1 list into a dataframe with category frequencies in each mode as separate cols
    vN <- as.matrix(vapply(vpt1, function(x) x$Freq, FUN.VALUE = numeric(length(nams))))
    rownames(vN) <- vpt1[[1]]$x
    colnames(vN) <- paste("Octave", fit$Data$octave, sep = "_")
  }
  
  if (plot_modes == T){
    if (!is.null(categ)){
      
      if (is.null(col.statu)) {cols = rainbow(length(unique(dat$status)))} else {
        cols = col.statu
      }
      yLim = max(colSums(vN)) * 1.10
      barplot(vN, col = cols, names.arg = fit$Data$octave, ylim = c(0, yLim))
      legend("topright", legend = rownames(vN), pch = 15, col = cols)
      points(fit, pch = 16, col = "black")
      
    } else {
      cols <- ifelse(fit$Data$octave %in% peaks, "red", "grey")
      yLim = max(predict(fit)) * 1.10
      barplot(fit, col = cols, ylim = c(0, yLim))
      points(fit, pch = 16, col = "black")
    }
  }
  
  if (is.null(categ)) {
    res <- list("Peak_locations" = peaks, "Species_per_octave" = HE)
  } else {
    res <- list("Peak_locations" = peaks, "Species_per_octave" = HE, 
                "Summary_table" = vN)
  }
  
  class(res) <- c("deconstruct")
  return(res)
}






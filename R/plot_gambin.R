#' @importFrom graphics barplot lines points
#' @export
barplot.gambin = function(height, ...) {
  
  ylim = max(predict(height), height$Data$species)
  defaults = list(xlab = "Octaves", ylab = "Number of species", 
                  ylim = c(0, 0.5 + ylim), names.arg = height$Data$octave)
  args = list(...)
  for(default in names(defaults)) {
    if(is.null(args[[default]])) args[[default]] = defaults[[default]]
  }
  do.call(barplot, c(list(height$Data$species), args))
  return(invisible(data.frame(x = height$Data$octave, y = height$Data$species)))  
}

#' @export
lines.gambin = function(x, ...) {
  midpoints = barplot(x$Data$species, plot = FALSE)
  lines(midpoints, predict(x), ...)
  return(invisible(data.frame(x = midpoints, y = predict(x))))
}

#' @export
points.gambin = function(x, ...) {
  midpoints = barplot(x$Data$species, plot = FALSE)
  points(midpoints, predict(x), ...)
  return(invisible(data.frame(x = midpoints, y = predict(x))))
}

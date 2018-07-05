#' @export
print.gambin = function(x, ...) 
{
  cat("\n")
  cat("Gambin distribution fit\n")
  cat("Data: ")
  cat(x$Dataname)
  cat("\n\n")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = 4), print.gap = 2L, quote = FALSE)
  cat("\n")
}

#' @export
print.summary.gambin = function(x, ...) 
{
  cat("Gambin distribution fit\n")
  cat("Data: ")
  cat(x$Dataname)
  cat("\n\n")
  
  alpha = round(x$alpha, 3)
  if(length(alpha) == 1L) {     
    vals <- data.frame(Estimated = alpha, CI95_low = x$ConfInt95[1], CI95_high = x$ConfInt95[2])
    rownames(vals) = paste0("alpha", seq_along(x$Alpha))
    cat("Coefficients:\n")
    print(vals)
    cat("\n")
  } else {
    if (length(x$ConfInt95) == 1){
      cat("No confidence intervals calculated", "\n\n")
    } else {
      cat("Confidence intervals:", "\n")
      print(x$ConfInt95)
      cat("\n\n")
    }
    cat("Alpha:\t")
    cat(x$alpha)
    cat("\n\n")
  }

  cat("MaxOctave:\t")
  cat(floor(x$octaves))
  
  cat("\n\n")
  cat(paste("Chi-squared fit:", 
            paste("X^2 =", round(x$ChiSq$statistic, 3)), 
            paste("df =", attr(x$logLik, "df")), 
            paste("p-value =", round(x$ChiSq$p.value, 3)), 
            "\n", sep = "\t" ))
  cat("\n")
}

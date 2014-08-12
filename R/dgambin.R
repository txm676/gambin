dgambin <-
function(alpha, maxoctave)
{
  # a function to give the probability distribution of the GamBin in 
  Pk <- function(alpha, k, octmax) # k is in 0:(nOct-1)
  {
    # calculates the 'fitness' distribution of species for a given alpha
    Gj <- function(alpha, j) # j is an integer in 1:100
    {
      qG99 = qgamma(0.99, alpha, 1) / 100
      b1 = (j - 1) * qG99
      b2 = j * qG99
      
      return((pgamma(b2, alpha, 1) - pgamma(b1, alpha, 1)) / 0.99)
    }
    
    s <- 0
    for (j in 1:100)
      s[j] <- choose(octmax, k) * (j / 100) ^ k * (1 - j/100) ^ (octmax - k) * Gj(alpha, j)
    return(sum(s))
  }
  
  ret <- sapply(0:maxoctave, Pk, alpha = alpha, octmax = maxoctave)
  names(ret) <- 0:maxoctave
  ret
}

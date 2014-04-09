#' Generate random number(s) from a beta distribution with specified mean and standard deviation
#' 
#' \code{rbeta2} converts the mean and stdev to appropriate shape parameters and calls \code{rbeta}.
#' It has similar functionality to \code{betaval} but is faster.
#' 
#' @param n number of random deviates desired
#' @param mean mean of the distribution (must be in (0,1))
#' @param sdev standard deviation of the distribution. Limits depend on the value of \code{mean}
#' 
#' @details All arguments must be scalars.
#' @return a vector of length \code{n} with the random numbers
#' @seealso \code{\link{rbeta}} and \code{\link{betaval}}
rbeta2 <- function(n, mean, sdev) {
    stopifnot(is.numeric(n), is.numeric(mean), is.numeric(sdev), 
              length(n) == 1, length(mean) ==1, length(sdev) == 1)
  
  if (mean > 1 || mean < 0) { # Ensure parameter is within bounds
    stop("Attempted to use a mean for the beta distribution that is outside of [0,1]")
  } else if (sdev == 0) { # No randomness - return n copies of mean
    warning("Requesting random variates with no variation. Returning replicates of mean")
    return(rep(mean, n))
  } else {
    sigma2 <- sdev^2
    if (sigma2 >= (1 - mean) * mean) { # Ensure parameter is within bounds
      stop("Standard deviation too high for beta distribution")
    }
    # Convert mean and variance to shape parameters
    temp <- ((mean * (1 - mean)/(sigma2)) - 1)
    a <- mean * temp
    b <- (1 - mean) * temp
    # Get the random variates
    return(rbeta(n, a, b))
  }
}
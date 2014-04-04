#' Generate random number(s) from a beta distribution with specified mean and standard deviation
#' 
#' \code{rbeta2} converts the mean and stdev to appropriate shape parameters and calls \code{rbeta}.
#' It has similar functionality to \code{betaval} but is faster.
#' 
#' @param n number of random deviates desired
#' @param mean mean of the distribution (must be in (0,1))
#' @param sdev standard deviation of the distribution. Limits depend on the value of \code{mean}
#' 
#' @return a vector of length \code{n} with the random numbers
#' @seealso \code{\link{rbeta}} and \code{\link{betaval}}
rbeta2 <- function(n, mean, sdev) {
  if (mean > 1 || mean < 0) { # Ensure parameter is within bounds
    stop("Please select a mean beta between 0 and 1")
  }
  if (sdev == 0) { # No randomness
    bb <- rep(mean, n) 
  } else {
    vari <- sdev^2
    if (vari >= (1 - mean) * mean) { # Ensure parameter is within bounds
      stop("Standard deviation too high for beta distribution")
    }
    # Convert mean and variance to shape parameters
    vv <- mean * ((mean * (1 - mean)/(vari)) - 1)
    ww <- (1 - mean) * ((mean * (1 - mean)/(vari)) - 1)
    # Get the random variates
    bb <- rbeta(n, vv, ww)
  }
}
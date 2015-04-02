#' Calculate a CDF of quasi-extinction risk for the Ricker model
#' @param r the low-density growth rate
#' @param K the equilibrium density
#' @param sigma2 the variance of the log growth rates
#' @param time_horizon number of years to simulate
#' @param N0 initial population size
#' @param nsims number of simulations
#' @param Nx quasi-extinction threshold
#' @param d the difference between \code{log(N0)} and \code{log(Nx)}. Either specify d or specify Nx
#'   
#' @return A data frame with columns \code{Time} and \code{CDF}
extRiskRicker <- function(r, K, sigma2, time_horizon, N0, nsims=1000, Nx=NULL, d=log(N0/Nx) ) {

  #' Some error checking:
  if (length(d)==0) {
    stop("Either specify Nx or specify d")
  }
  if (is.null(Nx)) {
    Nx <- N0/exp(d)
  }
  if (Nx >= N0) {
    stop("Nx must be less than N0 (or d must be positive)")
  }
  replicated_simulations <- simulateRicker(r, K, sigma2, N0, time_horizon, nsims)
  quasi_extinctions <- replicated_simulations < Nx
  cumulative_quasi_extinctions <- apply(quasi_extinctions, 2, cummax)
  quasi_extinction_CDF <- apply(cumulative_quasi_extinctions, 1, mean)
  return(data.frame(Time=1:time_horizon, CDF=quasi_extinction_CDF))
}
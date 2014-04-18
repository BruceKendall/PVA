#' Calculate a CDF of quasi-extinction risk for the SEG model
#' @param mu the mean log growth rate
#' @param sigma2 the variance of the log growth rates
#' @param time_horizon number of years to simulate
#' @param nsim number of replicate simulations for estimating the CDF
#' @param N0 initial population size
#' @param Nx quasi-extinction threshold
#' @param d the difference between \code(log(N0)) and \code(log(Nx)). Either specify d or specify Nx
#'   and N0
#'   
#' @return a data frame with columns \code{Time} and \code{CDF}
extRiskSEG <- function(mu, sigma2, time_horizon, nsims=1000, N0=NULL, Nx=NULL, d=log(Nx/N0) ) {
  if (is.null(d)) {
    stop("Either specify N0 and Nx or specify d")
  }
  if (is.null(N0)) {
    N0 <- 1
  }
  if (is.null(Nx)) {
    Nx <- N0/exp(d)
  }
  if (Nx >= N0) {
    stop("Nx must be less than N0 (or d must be positive)")
  }
  replicated_SEG_simulations <- replicate(nsims, simulateSEG(mu, sigma2, N0, time_horizon))
  quasi_extinctions <- replicated_SEG_simulations < Nx
  cumulative_quasi_extinctions <- apply(quasi_extinctions, 2, cummax)
  quasi_extinction_CDF <- apply(cumulative_quasi_extinctions, 1, mean)
  return(data.frame(Time=1:time_horizon, CDF=quasi_extinction_CDF))
}
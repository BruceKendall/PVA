#' Simulate a time series from the stochastic exponential growth model
#' @param mu the mean log growth rate
#' @param sigma2 the variance of the log growth rates
#' @param N0 initial population size
#' @param time_horizon number of years to simulate
#' 
#' @return vector of length \code{time_horizon} with the simulated abundances
simulateSEG <- function (mu, sigma2, N0, time_horizon) {
  sigma <- sqrt(sigma2)
  Nt <- numeric(time_horizon+1)
  Nt[1] <- N0
  lambdat <- exp(rnorm(time_horizon, mu, sigma))
  for (Time in 1:time_horizon) {
    Nt[Time+1] <- Nt[Time] * lambdat[Time]
  }
  Nt <- Nt[-1]
  Nt
}

#' Simulate one or more time series from the stochastic Ricker model
#' @param r the low-density growth rate
#' @param K the equilibrium density
#' @param sigma2 the variance of the log growth rates
#' @param N0 initial population size
#' @param time_horizon number of years to simulate
#' @param nsims number of simulations
#' @param keep_N0 (logical) whether to return \code{N0} at the beginning of the time series
#' 
#' @return If \code{nsim=1} (the default), a vector of length \code{time_horizon} (or \code{time_horizon+1} if \code{keep_N0=TRUE}) with the simulated abundances. Otherwise, a matrix with \code{time_horizon} (or \code{time_horizon+1} if \code{keep_N0=TRUE}) rows and \code{nsim} columns.
#' 
#' @details This function simulates the Ricker discrete time model with environmental stochasticity, defined as follows
#' \deqn{N(t+1) = N(t) exp[r (1 - N(t)/K) + \epsilon(t)]}, 
#' where
#' \deqn{\epsilon(t) ~ Normal(0, \sigma)}
#' 
#' @examples
#' simulateRicker(2, 100, 0.01, 20, 10)
#' simulateRicker(2, 100, 0.01, 20, 10, 5, keep_N0=TRUE)
simulateRicker <- function (r, K, sigma2, N0, time_horizon, nsims=1, keep_N0=FALSE) {
  sigma <- sqrt(sigma2)
  # We keep track of both abundance and log abundance to reduce the number
  #   of transendental function evaluations
  Nt <- matrix(0, time_horizon+1, nsims)
  Xt <- matrix(0, time_horizon+1, nsims)
  Nt[1,] <- N0
  Xt[1,] <- log(N0)
  
  # Pre-calculate the epsilons for efficiency
  epsilont <- matrix(rnorm(time_horizon*nsims, 0, sigma), time_horizon, nsims)
  
  for (Time in 1:time_horizon) {
    Xt[Time+1,] <- Xt[Time,] + r * (1 - Nt[Time,]/K) + epsilont[Time,]
    Nt[Time+1,] <- exp(Xt[Time+1,])
  }
  
  if (!keep_N0) { # Drop the first row
    Nt <- Nt[-1,]
  }
  
  if (nsims==1) { # Return a vector
    as.vector(Nt)
  } else { # Return the matrix
    Nt
  }
}

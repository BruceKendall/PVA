\name{simulateSEG}
\alias{simulateSEG}
\title{Simulate a time series from the stochastic exponential growth model}
\usage{
simulateSEG(mu, sigma2, N0, time_horizon)
}
\arguments{
  \item{mu}{the mean log growth rate}

  \item{sigma2}{the variance of the log growth rates}

  \item{N0}{initial population size}

  \item{time_horizon}{number of years to simulate}
}
\value{
vector of length \code{time_horizon} with the simulated
abundances
}
\description{
Simulate a time series from the stochastic exponential
growth model
}


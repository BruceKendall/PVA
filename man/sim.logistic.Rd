\name{sim.logistic}

\alias{sim.logistic}

\title{Simulate the Stochastic Logistic Model}

\description{Simulate the Stochastic logistic Model}

\usage{
sim.Ricker(N0, param, Tmax, nsim, method=c("normal","resample"), resids=NULL)
}

\arguments{
\item{N0}{Starting population size}
\item{param}{Vector of parameters for the model.  First element is R, second is
    K, and the third (used only if \code{method="normal"}) is the variance of
    the epsilons.}
\item{Tmax}{Number of years to simulate}
\item{nsim}{Number of replicate simulations}
\item{method}{Treat epsilon as normally distributed random numbers (\code{"normal"})
    or resample from the observed residuals (\code{"resample"}).  The default is
    \code{"normal"}}
\item{resids}{The residuals from which to resample to simulate epsilon.  
    Used only if \code{method="resample"}}
}

\details{
This is pretty straightforward.
}

\value{
A matrix with \code{Tmax} rows and \code{nsim} columns.  Each column is the time
series from a replicate simulation.
}

\references{
    Morris, W. F. and Doak D. F. (2002) \emph{Quantitative Conservation Biology:
        Theory and Practice of Population Viability Analysis}.
        Sunderland: Sinauer Associates.
}


\author{Bruce E. Kendall (\email{kendall@bren.ucsb.edu})}


\examples{
## Generate 10 replicate simulations 50 years long
sim.Ricker(100, param=c(2,120,.2), Tmax=50, nsim=10)
}

\keyword{models}

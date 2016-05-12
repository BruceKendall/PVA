\name{sim.delay}

\alias{sim.delay}

\title{Simulate the Stochastic Delay-Ricker Model}

\description{Simulate the Stochastic delay-Ricker Model}

\usage{
sim.delay(N0, param, Tmax, nsim, method=c("normal","resample"), resids=NULL)
}

\arguments{
\item{N0}{Starting population size. A two-element vector, containing population
    size at times 0 and -1}
\item{param}{Vector of parameters for the model.  First element is a0, second is
    a1, third is a2, 
    and the fourth (used only if \code{method="normal"}) is the variance of
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
sim.delay(c(130,100), param=c(0.5,-0.001,-0.0001,.2), Tmax=50, nsim=10)
}

\keyword{models}

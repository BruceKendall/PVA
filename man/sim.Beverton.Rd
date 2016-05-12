\name{sim.Beverton}

\alias{sim.Beverton}

\title{Simulate the Stochastic Beverton-Holt Model}

\description{Simulate the Stochastic Beverton-Holt Model}

\usage{
sim.Beverton(N0, param, Tmax, nsim, method=c("normal","resample","DS"), resids=NULL)
}

\arguments{
\item{N0}{Starting population size}
\item{param}{Vector of parameters for the model.  First element is R, second is
    a, and the third (used only if \code{method="normal"}) is the variance of
    the epsilons.  If \code{method="DS"}, then the third parameter is the intercept
    and the fourth is the slope of the regression of squared residuals on 1/N.}
\item{Tmax}{Number of years to simulate}
\item{nsim}{Number of replicate simulations}
\item{method}{Treat epsilon as normally distributed random numbers (\code{"normal"}),
    resample from the observed residuals (\code{"resample"}), or have the variance be
    inversely related to population size to simulate demographic stochasticity 
    (\code{method="DS"}).  The default is \code{"normal"}.}
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
sim.Beverton(100, param=c(2,0.002,.1), Tmax=50, nsim=10)
}

\keyword{models}

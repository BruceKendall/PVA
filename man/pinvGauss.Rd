\name{pinvGauss}

\alias{pinvGauss}

\title{Inverse Gaussian distribution}

\description{
    Cumulative distribution function for the inverse Gaussian distribution, for
        use with count-based density-independent PVA.
}

\usage{
    pinvGauss(T, mu, s2, d)
}

\arguments{
    \item{T}{Time point(s) at which to evaluate the distribution.  May be a vector.}
    \item{mu}{Mean growth rate of the population.}
    \item{s2}{Variance of the population growth rate.}
    \item{d}{Natural log of the ratio of starting population size to quasi-extinction
            population size.  In the language of Morris and Doak (2002), 
            d = log(Nc) - log(Nx).}
}

\details{
Calculates the cumulative distribution function for the inverse Gaussian 
distribution, using the parameterization of Eq. (3.5) in Morris and Doak 
(2002).  It is equivalent to the Matlab program \code{extcdf.m} on p. 80
of that book.  This can be interpreted as the probability of quasi-extinction
by a given time for a population exhibiting density-independent growth with 
the mean and variance of the instantaneous growth rate, r = log(lambda), 
given by mu and s2.  See Chapt. 3 of Morris and Doak for details.
}

\value{
A vector of the same length as T, containing the CDF of the inverse Gaussian 
distribution at each value of T.
}

\references{
    Morris, W. F. and Doak D. F. (2002) \emph{Quantitative Conservation Biology:
        Theory and Practice of Population Viability Analysis}.
        Sunderland: Sinauer Associates.
}

\note{
    This differs from the inverse Gaussian distributions in the \pkg{statmod} and
    \pkg{SuppDists} libraries, not only in the parameterization, but also in that it 
    allows evaluation of the incomplete inverse Gaussian distribution (the
    integral of the pdf is less than 1) that occurs when mu > 0.
}

\author{Bruce E. Kendall (\email{kendall@bren.ucsb.edu})}

\seealso{
    \code{\link{count.DI.PVA}}
}

\examples{
## Load the grizzly bear data
data(grizzly)

## Calculate the log-growth rates, and calculate their mean and variance
r <- log(grizzly$N[-1]/grizzly$N[-39])
mu <- mean(r); s2 <- var(r)

## Calculate the cumulative quasi-extinction risk (with threshold at 20) at
##  years 1-100
pinvGauss(1:100, mu, s2, log(grizzly$N[39]/20))
}

\keyword{distribution}

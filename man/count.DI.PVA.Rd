\name{count.DI.PVA}

\alias{count.DI.PVA}

\title{Count-Based Density-Independent PVA}

\description{
Implements the method of Dennis et al. (1991), as described in Morris and Doak
(2002)
}

\usage{
    count.DI.PVA(Nt, Yeart=seq(along=Nt), Nc=Nt[length(Nt)], Nx=1, year.max=100, 
        nboot=0, alpha=0.05, meas.err=c("none","simple"), Nt.se=NULL, ...) 
}

\arguments{
    \item{Nt}{Vector containing the time series of population abundance estimates}
    \item{Yeart}{Vector of census times; must be the same length as \code{Nt}.  Defaults to
        a sequence of integers}
    \item{Nc}{Current population size.  Defaults to the last value of \code{Nt}}
    \item{Nx}{The quasi-extinction threshold.  Defaults to 1}
    \item{year.max}{The number of years over which to do the projection.  Defaults
                to 100.}
    \item{nboot}{Number of bootstrap replicates to use for putting confidence
        intervals on the extinction risk calculations.  Defaults to 0.}
    \item{alpha}{One minus the desired confidence limits.}
    \item{meas.err}{Is there measurement error in the census data?  If not, set 
        \code{meas.err="none"} (the default); \code{meas.err="simple"} uses eq. (5.5)
        in Morris and Doak (2002) to correct the estimatee of the sampling variance.}
    \item{Nt.se}{Standard errors of each of the population estimates. Must be specified
        if \code{meas.err="simple"}.}
    \item{...}{Additional plotting parameters.}
}

\details{
Implements the method of Dennis et al. (1991), as described in Morris and Doak
(2002).  If the censuses are not evenly spaced in time, then the parameters
(mu and sigma2) are estimated using the regression method.

If \code{nboot=0}, then all the components related to the confidence intervals are set
to \code{NULL}.  Otherwise, nonparametric bootstrap is performed: either the rt values 
or the (x,y) pairs for the regression are sampled with replacement for each 
bootstrap replicate.

The confidence intervals are pointwise estimates at each projection year, and so
might not represent any actual extinction curve.

When correcting for sampling errors, the smallest allowable value of sigma2 is 0.001.
Smaller values result in numerical errors in the inverse Gaussian function.  This correction 
only works if you have a continuous time series -- I'm not sure how to do the correction
when using the regression approach to estimate mu and sigma2.
}

\value{
    Returns an object of class \code{\link{ext.risk}}.  Also displays a plot of the 
    CDF of extinction risk; if \code{nboot>0} then the plot also shows the confidence 
    intervals.
}

\references{
    Dennis, B., Munholland, P. L. and Scott, J. M. (1991) Estimation of growth
        and extinction parameters for endangered species.  \emph{Ecological
        Monographs} \bold{61}, 115--143.
    
    Morris, W. F. and Doak D. F. (2002) \emph{Quantitative Conservation Biology:
        Theory and Practice of Population Viability Analysis}.
        Sunderland: Sinauer Associates.
}

\note{
    This function assumes that the units of time are years, and that there is at
    most one census per year.  If these are violated, I make no guarantee as to the 
    meaningfulness of the results....
    
    There is currently no testing of the sensibility of the input parameters, so
    if you get an arcane error message check that your inputs make sense and have
    correct dimensions.
}

\author{Bruce E. Kendall (\email{kendall@bren.ucsb.edu})}

\seealso{
    \code{\link{pinvGauss}} for the important part of the computation; 
    for a description of the components of the 
    returned object see \code{\link{ext.risk}}.
}

\examples{
## Load the grizzly bear data
data(grizzly)

## Run the PVA with the defaults, plotting risk on a log scale, and print the result
count.DI.PVA(grizzly$N, log='y')

## Specify some non-default parameters, and get confidence intervals
grizzly.ext <- count.DI.PVA(grizzly$N, Nx=20, Nc=50, nboot=1000, year.max=500)
plot(grizzly.ext)
grizzly.ext
}

\keyword{models}

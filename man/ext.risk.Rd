\name{ext.risk}

\alias{ext.risk}

\title{Create an ext.risk Object}

\description{
Create an object of class \code{ext.risk}.
}

\usage{
ext.risk(extrisk, year, nboot, ext.boot, alpha, ext.ci, model.type, model.params)
}

\arguments{
\item{extrisk}{The cumulative probability of quasi-extinction projected at various
    times in the future}
\item{year}{The times associated with the risks in \code{extrisk}}
\item{nboot}{Number of bootstrap replicates used to calculate the confidence intervals}
\item{ext.boot}{Extinction risk curves for each of the bootstrap replicates}
\item{alpha}{One minus the confidence level}
\item{ext.ci}{Confidence intervals of cumulative extinction risk}
\item{model.type}{Name of the model used to perform the PVA}
\item{model.params}{Estimated parameters of the model used to perform the PVA}
}

\details{
If \code{extrisk} and \code{year} are vectors of length n, then \code{ext.boot} has dimension 
\code{(n, nboot)} and \code{ext.ci} has dimension \code{(n, 2)}.  Each column
of \code{ext.boot} represents a single bootstrap replicate.

If there are no bootstrap replicates, set \code{nboot, ext.boot, alpha}, and 
\code{ext.ci} all to \code{NULL}.
}

\value{
An object of class ext.risk
}


\note{
This function was developed in order to document the \code{ext.risk} class.  In the
normal course of events, there should be no need for the user to call this 
function.
}

\author{Bruce E. Kendall (\email{kendall@bren.ucsb.edu})}

\seealso{
\code{\link{print.ext.risk}} and \code{\link{plot.ext.risk}} for printing and
plotting \code{ext.risk} objects; these objects are generated, for example, 
by \code{\link{count.DI.PVA}}
}

\examples{
## Load the grizzly bear data and run a DI PVA
data(grizzly)
grizzly.ext <- count.DI.PVA(grizzly$N, Nx=20, Nc=50, nboot=1000)

## Print the confidence intervals for all 100 years
cbind(grizzly.ext$year, grizzly.ext$ext.ci)

## Plot a selection of bootstrap replicates
matplot(grizzly.ext$year, grizzly.ext$ext.boot[,1:100], type='l', col='black', 
    lty=1, xlab="Year", ylab="Cumulative extinction probability")
}

\keyword{classes}

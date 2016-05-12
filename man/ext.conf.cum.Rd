\name{ext.conf.cum}

\alias{ext.conf.cum}

\title{Extract Confidence Levels for Extinction Risk}

\description{
For a given projection time, extract the cumulative confidence levels for 
the extinction risk (confidence that the true extinction risk is at or below
a certain level).
}

\usage{
ext.conf.cum(x, year=length(x$year)) 
}

\arguments{
\item{x}{An ext.risk object.  \code{x$ext.boot} must not be \code{NULL}}
\item{year}{The projection year at which to evaluate the confidence levels.
    Defaults to the last year of the PVA}
}

\details{
}

\value{
A data frame with two columns:


\emph{Extinction.risk}: The extinction risk

\emph{Cumulative.confidence}: The confidence that the extinction risk is at or
    below the given value
}


\author{Bruce E. Kendall (\email{kendall@bren.ucsb.edu})}

\examples{
## Do a PVA of the checkerspot data
data(checkerspot)
ricker.ext <- count.DD.PVA(checkerspot$N, model="Ricker", nboot=100)

## Plot the cumulative confidence curve for the 20-year projection
plot(ext.conf.cum(ricker.ext, 20))
}

\keyword{models}

\name{AICc}

\alias{AICc}

\title{Corrected AIC}

\description{
Calculates the corrected (for small sample size) Akaike Information Criterion
for \code{lm} and \code{nls} objects.
}

\usage{
AICc(x)
}

\arguments{
\item{x}{An \code{lm} or \code{nls} object}
}

\details{
See eq. (4.8) in Morris and Doak (2002).
}

\value{
The corrected AIC of the model.
}

\references{
    Morris, W. F. and Doak D. F. (2002) \emph{Quantitative Conservation Biology:
        Theory and Practice of Population Viability Analysis}.
        Sunderland: Sinauer Associates.
}


\author{Bruce E. Kendall (\email{kendall@bren.ucsb.edu})}


\examples{
## Load the Jasper Ridge checkerspot butterfly data
data(checkerspot)
rt <- diff(log(checkerspot$N))
Nt <- checkerspot$N[-dim(checkerspot)[1]]

## Fit a Ricker model and extract the corrected AIC
ricker.nls <- nls(rt ~ mu*(1-Nt/K), start=list(mu=1,K=500))
AICc(ricker.nls)
}

\keyword{models}

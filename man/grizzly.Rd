\name{grizzly}

\docType{data}

\alias{grizzly}

\title{Yellowstone Grizzly Bear Population Data}

\description{
Estimated number of adult female grizzly bears in the Greater Yellowstone
population, each year from 1959 to 1997
}

\usage{
data(grizzly)
}

\details{
Quoting Morris and Doak (2002), p. 65:

This population is currently designated as threatened by the United States Fish
and Wildlife Service, and is completely isolated from other grizzly bear
populations.  Each year, bear biologists count the number of unique female bears
with cubs (offspring in the first year of life) in the entire Yellowstone population.
Because the litters of one to three cubs remain closely associated with their 
mothers, females with cubs are the most conspicuous, and therefor most reliably 
censused, element of the population.  Censuses were originally performed by 
observing bears at park garbage dumps, but they havbe been conducted by aerial
survey since the closure of the dumps in 1970--1971.  The counts of females with
cubs are used to estimate the total number of adult females in the population.  
Specifically, the number of adult females in year \eqn{t} is estimated as the 
sum of the observed numbers of females with cubs in years \eqn{t}, \eqn{t+1}, and
\eqn{t+2} (Eberhardt et al. 1986).  The logic underlying this estimate is that the
interval between litters produced by the same mother is at least three years, so
that females with cubs observed in years \eqn{t+1} and \eqn{t+2} could not have 
been the same individuals that were observed with cubs in year \eqn{t}.  Yet, if
they were observed in later years, they must have been alive in year \eqn{t} 
(although they may not have been \emph{adults} in year \eqn{t}, which introduces
some error).
}

\format{
A data frame with three columns:
\describe{
\item{i}{Census index}
\item{yr}{Census year}
\item{N}{Estimated number of bears}
}
}

\source{
Morris and Doak (2002), table 3.1.  Data originally from Eberhardt et al.
(1986) and Haroldson (1999).
}

\references{
Eberhardt, L. L., Knight, R. R. and Blanchard, B. (1986) Monitoring grizzly bear
population trends.  \emph{Journal of Wildlife Management} \bold{50}, 613--618.

Haroldson, M. (1999) Bear monitoring and population trend: unduplicated females.
Pages 3--17 in Schwartz, C. C. and Haroldson, M. A. (eds.) \emph{Yellowstone
Grizzly Bear Investigations: Annual Report of the Interagency Grizzly Bear Study
Team, 1999}. US Geological Survey (http://www.nrmsc.usgs.gov/research/igbst-home.htm).

    Morris, W. F. and Doak D. F. (2002) \emph{Quantitative Conservation Biology:
        Theory and Practice of Population Viability Analysis}.
        Sunderland: Sinauer Associates.
}

\author{Bruce E. Kendall (\email{kendall@bren.ucsb.edu})}

\keyword{datasets}

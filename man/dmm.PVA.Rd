\name{dmm.PVA}

\alias{dmm.PVA}

\title{Run a Demographic PVA}

\description{
Perform a demographic PVA, incorporating environmental stochasticity, 
demographic stochasticity, or both.
}

\usage{
dmm.PVA(x, Nc, Nx, year.max=100, nsim=1000, 
    ES=c("MatrixDraw","ParDraw","none"), DS=FALSE, Varsamp=c("none","White"), 
    fecmod=c("logN","sBeta"), fecmax=1, ...) }

\arguments{
\item{x}{A demographic matrix model, created with \code{dmm}}
\item{Nc}{Current population size.  Must be a vector with a value for the number of 
    individuals in each class}
\item{Nx}{The quasi-extinction threshold.  A scalar; extinction is considered to
    have occured if the total population size falls below this value.}
\item{year.max}{The number of years over which to do the projection.  Defaults
                to 100.}
\item{nsim}{Number of replicate simulations.  Defaults to 1000}
\item{ES}{Type of environmental stochasticity model.  Defaults to "MatrixDraw". 
See "Details" for more information.}
\item{DS}{Logical flag for whether to include demographic stochasticity.}
\item{Varsamp}{Sampling correction for the variance estimates.  Defaults to "none".
See "Details" for more information.}
\item{fecmod}{Distribution to use for environmental stochasticity in fertility:
log-normal ("logN") or stretched Beta ("sBeta").  Defaults to "logN".}
\item{fecmax}{Maximum possibility, for use when fecmod="sBeta".}
}

\details{
When \code{ES = "MatrixDraw"}, the annual matrices are drawn at random for use in 
the simulation.  When \code{ES = "ParDraw"}, the individual demographic rates
are drawn at random, based on the variance-covariance structure in the annual 
matrices.  The construction of correlated random variables from arbitrary 
distributions follows the approach described in Chapter 8 of Morris and Doak
(2002).  Annual variation in survival and growth rates is treated as beta-distributed;
annual variation in fertility is either log-normal or stretched-beta, depending
on the value of \code{fecmod}.  Animals with at most one offspring in a given
year, such as most seabirds and many large mammals, should have \code{fecmod = 
"sBeta"} and \code{fecmax=1}.

When \code{DS = FALSE} (the default), there is no demographic stochasticity and 
the random matrix is simply multiplied by the existing population vector to get
the simulated population vector in the following year.  When \code{DS = TRUE} then 
the number of individuals surviving and growing are drawn from binomial distributions,
and the number of newborns are drawn from Poisson distributions.

When estimating the among year variances in the parameters, you can use
White's correction for sampling error in the variance estimate
(White 2000; see Morris and Doak 2002 Ch. 8) by setting \code{Varsamp = "White"}.
}

\value{
An object of class ext.risk
}

\references{
   Morris, W. F. and Doak D. F. (2002) \emph{Quantitative Conservation Biology:
        Theory and Practice of Population Viability Analysis}.
        Sunderland: Sinauer Associates.
        
Whte, G. C. (2000) Population viability analysis: data requirements and essential
    analyses.  Pp. 288--331 in Boitani, L. and Fuller, T. K. (eds.) \emph{Research
    Techniques in Animal Ecology: Controversies and Consequences}.
    New York: Columbia University Press. 
}

\note{
This runs relatively slowly, so during model development you will want to set
\code{year.max} and \code{nsim} below their default values.

In general, this code will \emph{not} work correctly for a model with a seedbank
and a pre-breeding census (i.e., with two rows of reproduction transitions.  Until
this is fixed, use a post-breeding census if you have a seed bank.

Provision for bootstrapping across years will be added in the future, as will 
additional options for demographic stochasticity in fertility.
}

\author{Bruce E. Kendall (\email{kendall@bren.ucsb.edu})}

\seealso{
\code{\link{dmm}}  for creating \code{dmm} objects; \code{\link{dmm.iter}} and 
\code{\link{dmm.asymp}} for deterministic analysis of matrix models.
}

\examples{
## Read in the Cereus data
data(Cereus.Grow)
data(Cereus.Fruit)

## Extract annual demography from the cereus data
cereus2.dmm <- dmm(Cereus.Grow, fec=Cereus.Fruit, p0=0.03448, type="size",
    census="pre", data.type="aggregated",
    classes=c("<1","1-4","4-8","8-16","16-32",">32"),annual=T,years=1988:1992)

# Estimated extinction risk, resampling matrices to simulate Environmental
#   stochasticity
cereus1.ext <- dmm.PVA(cereus2.dmm, Nc=c(10,11,19,16,25,15), Nx=50, 
    year.max=100, nsim=1000, ES="MatrixDraw", DS=FALSE)

# Estimated extinction risk, drawing random parameters to simulate Environmental
#   stochasticity, using the White correction for sampling uncertainty in the
#   variance, and incorporating demographic stochasticity
cereus2.ext <- dmm.PVA(cereus2.dmm, Nc=c(10,11,19,16,25,15), Nx=50, 
    year.max=100, nsim=1000, ES="ParDraw", Varsamp="White",
    DS=TRUE, fecmod="logN")

# Estimated extinction risk, with demographic stochasticity alone
cereus3.ext <- dmm.PVA(cereus2.dmm, Nc=c(10,11,19,16,25,15), Nx=50, 
    year.max=100, nsim=1000, ES="none", DS=TRUE)
}

\keyword{classes}

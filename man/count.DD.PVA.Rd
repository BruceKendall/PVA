\name{count.DD.PVA}

\alias{count.DD.PVA}

\title{Count-Based Density-Dependent PVA}

\description{
Simulation-based method, as described in Morris and Doak
(2002)
}

\usage{
count.DD.PVA(Nt, Yeart=seq(along=Nt), Nc=Nt[length(Nt)], Nx=1, 
    year.max=100, nsim=1000, nboot=0, alpha=0.05, 
    boot.type=c("nonparametric","parameteric"),
    model=c("Ricker","theta","logistic","Beverton","delay"), 
    enoise=c("normal","resample","DS"),
    par.start=NULL, meas.err=FALSE, ...) 
}

\arguments{
    \item{Nt}{Vector containing the time series of population abundance estimates.
        \emph{Must be yearly estimates.}}
    \item{Yeart}{Vector of census times; must be the same length as \code{Nt}.  Defaults to
        a sequence of integers}
    \item{Nc}{Current population size.  Defaults to the last value of \code{Nt}}
    \item{Nx}{The quasi-extinction threshold.  Defaults to 1}
    \item{year.max}{The number of years over which to do the projection.  Defaults
                to 100.}
    \item{nsim}{Number of replicate stochastic simulations}
    \item{nboot}{Number of bootstrap replicates to use for putting confidence
        intervals on the extinction risk calculations.  Defaults to 0.}
    \item{alpha}{One minus the desired confidence limits.}
    \item{boot.type}{Use parametric or nonparametric bootstrap?}
    \item{model}{Type of density-dependence model to use}
    \item{enoise}{In stochastic simulations, are the epsilons (environmenatal
        stochasticity) normally distributed (\code{enoise="normal"}), incorporate
        demograhic stochasticity (\code{enoise="DS"}) or resampled from the observed 
        residuals (\code{enoise="resample"})?}
    \item{par.start}{Initial estimates of the parameters, useful for difficult
        nonlinear least squares problems}
    \item{meas.err}{Not yet used.}
    \item{...}{Additional parameters to pass to plot or nls.}
}

\details{
Implements the simulation-based approach to count-based PVA, as described in 
Morris and Doak
(2002: ch. 4).  The censuses must be annual.  Currently implemented models are
Density-independent, Ricker, theta-logistic, logistic, Beverton-Holt, and delay-Ricker.

If \code{nboot=0}, then all the components related to the confidence intervals are set
to \code{NULL}.  Otherwise, nonparametric or parametric bootstrap is performed. For
nonparametric bootstrap (the default) the (rt, Nt) pairs are sampled with replacement
and the regression re-run to get the new parameter estimates.  For parametric
bootstrap, the parameters are sampled from the multivariate normal distribution
defined by the covariance matrix that is returned by lm or nls.

Both approaches have problems.  For models that are fit by nls, it is well-nigh
impossible to guarantee convergence of the fitting algorthithm on every 
nonparametric bootstrap replicate; this failure stops the whole program.  You can 
attempt to alleviate this by passing in control parameters (see \code{\link{nls}} and 
\code{\link{nls.control}}), but success is unlikely, especially for models like
the theta-logistic that are intrinsically ill-behaved.  See \bold{Examples} for a 
labour-intensive workaround.

In contrast, parametric bootstrap may generate biologically implausible parameter
values, such as negative values of K.  This creates programming as well as 
biological problems, for the resulting simulations tend to blow up.  Thus, the 
program discards any bootstrap replicates that have K < 0; this has the unfortunate
side effect that the mean of the bootstrapped parameter values does not equal
the original parameter estimates.

The confidence intervals are pointwise estimates at each projection year, and so
might not represent any actual extinction curve.

Demographic stochasticity is simulated by regressing the squared residuals on the inverse
of population size, and then using this relationship to generate a variance in the 
epsilons that depends inversly on population size.
}

\value{
    Returns an object of class \code{\link{ext.risk}}.  Also displays a plot of the 
    CDF of extinction risk; if \code{nboot>0} then the plot also shows the confidence 
    intervals.
}

\references{
    Morris, W. F. and Doak D. F. (2002) \emph{Quantitative Conservation Biology:
        Theory and Practice of Population Viability Analysis}.
        Sunderland: Sinauer Associates.
}

\note{
    This function assumes that the units of time are years, and that there is at
    most one census per year.  If these are violated, I make no guarantee as to the 
    meaningfulness of the results\ldots.
    
    There is currently no testing of the sensibility of the input parameters, so
    if you get an arcane error message check that your inputs make sense and have
    correct dimensions.
}

\author{Bruce E. Kendall (\email{kendall@bren.ucsb.edu})}

\seealso{
    \code{\link{sim.Ricker}}, \code{\link{sim.theta}}, and 
    \code{\link{sim.logistic}}, for the simulation models; 
    for a description of the components of the 
    returned object see \code{\link{ext.risk}}.
}

\examples{
## Load the checkerspot butterfly data
data(checkerspot)

## Run a PVA with the Ricker model
jrh.ext <- count.DD.PVA(checkerspot$N, Nx=20, nboot=100)

## Fit the theta-logistic model to the data
rt <- diff(log(checkerspot$N))
Nt <- checkerspot$N[-dim(checkerspot)[1]]
thetalog.nls <- nls(rt ~ mu*(1-(Nt/K)^theta), start=list(mu=1,K=500,theta=1))

## Run the corresponding PVA, using parametric bootstrap
jrh2.ext <- count.DD.PVA(checkerspot$N, Nx=20, nboot=100, boot.type="parametric",
    model="theta", par.start=coef(thetalog.nls))

## To generate nonparametric bootstrap reps for the theta logistic model,
##   first run the following two lines...
base.ext <- count.DD.PVA(checkerspot$N,Nx=20,nboot=0,model="theta")
boot.out <- as.matrix(base.ext$extrisk)  

## Then cut and paste the following four lines repeatedly until you have
##  100 reps.  Each column of boot.out will be a replicate
while (dim(boot.out)[2] < 100) {
    boot.ext <- count.DD.PVA(checkerspot$N,Nx=20,nboot=1,model="theta")
    boot.out <- cbind(boot.out, boot.ext$ext.boot)
}
 
}

\keyword{models}

\name{dmm}

\alias{dmm}

\title{Create a Demographic Matrix Model}

\description{
Create a demographic matrix model object
}

\usage{
dmm(x, fec=NULL, grow=NULL, p0=NULL, classes=NULL,
    type=c("age","stage","size","sfa"), census=c("pre","post"),
    terminal=FALSE, data.type=c("rates","aggregated","raw"), 
    fec.sampling=c("Poisson","binomial"), seedbank=FALSE,
    seed.germ=NULL, seed.surv=NULL, annual=F, years=NULL)

dmm(x, classes=NULL, census=c("pre","post"), type="age", terminal=F, annual=F, years=NULL)

dmm(x, grow, p0=NULL, classes, census=c("pre","post"), type, 
    data.type="rates", seedbank=FALSE, seed.germ=NULL, seed.surv=NULL)

dmm(x, fec, p0=NULL, classes=NULL, census=c("pre","post"), type,
    data.type, fec.sampling=c("Poisson","binomial"), seedbank=FALSE,
    seed.germ=NULL, seed.surv=NULL, annual=F, years=NULL)
}

\arguments{
These depend on the model being built; see "details".
}

\details{
In the first case, create an age-structured model:
\describe{
\item{\code{x}}{A numeric matrix containing the life table.  It contains three columns: the first
    is the age, the second is the survival of that age class, and the third is
    the fertility of that age class.  Must include information on age zero 
    individuals.  If \code{annual=TRUE}, then the first column must be the year, with the
    other columms shifted over to the right.}
\item{\code{classes}}{A character vector with the name of each age class.  If omitted,
    the class names are taken from the first column of \code{x}.}
\item{\code{census}}{Whether to construct a model that represents a pre-breeding census
    or a post-breeding census.  Defaults to \code{"pre"}.}
\item{\code{type}}{Indicates that we are doing an age-structured model.  May be omitted, 
    as this is the default.}
\item{\code{terminal}}{Indicates whether all individuals in the final age class die at the
    end of the year (\code{terminal = TRUE}), or whether they can remain in that 
    class indefinitely (\code{terminal = FALSE}, the default).}
\item{\code{annual}}{Logical.  Whether to extract annual models as well as the mean model.
    If set to \code{TRUE}, then \code{years} must be specified.}
\item{\code{years}}{A numeric list of the years for which annual models should be built.
    These must match values in the first column of \code{x}.}
}
    
In the second case, create a size-structured model using rate information:
\describe{
\item{\code{x}}{A numeric matrix containing the life table.  It contains three columns: 
    the first is an index of the size class, the second is the survival of that 
    size class, and the third is the fertility of that size class.  Must include 
    information on age zero  individuals \emph{only} if \code{census = "post"}. 
    If \code{annual=TRUE}, then the first column must be the year, with the
    other columms shifted over to the right.}
\item{\code{grow}}{A matrix of the growth transitions.  Element \{i,j\} of the matrix
    indicates the probability of growing from class j to class i, conditional
    on survival.  If \code{annual=TRUE}, then this must be a three-dimensional array,
    with \code{grow[,,i]} being the growth matrix for the i'th year.}
\item{\code{p0}}{Survival of newborns.  Must be specified if \code{census = "pre"}.}
\item{\code{classes}}{A character vector with the name of each age class.  Must be 
    specified.}
\item{\code{census}}{Whether to construct a model that represents a pre-breeding census
    or a post-breeding census.  Defaults to \code{"pre"}.}
\item{\code{type}}{Type of matrix.  Must be set to \code{"size"}.}
\item{\code{data.type}}{Specifies that we are providing rate information.  May be omitted.}
\item{\code{seedbank}}{Logical.  Indicates whether there is a persistent seed bank.  If 
    \code{TRUE}, then values must be specified for \code{seed.germ} and 
    \code{seed.surv}.  Defaults to \code{FALSE}.}
\item{\code{seed.germ}}{The fraction of surviving seeds in the seedbank that germinate.
    May be omitted if \code{seedbank = FALSE}.}
\item{\code{seed.surv}}{The annual survival of seeds in the seedbank.
    May be omitted if \code{seedbank = FALSE}.}
\item{\code{annual}}{Logical.  Whether to extract annual models as well as the mean model.
    If set to \code{TRUE}, then \code{years} must be specified.}
\item{\code{years}}{A numeric list of the years for which annual models should be built.
    These must match values in the first column of \code{x}.}
}

In the third case, create a size structured model using "aggregated" data:
\describe{
\item{\code{x}}{A matrix specifing the growth and survival data.  It has four columns:
    (1) the year; (2) an index for the size class at the start of the year; 
    (3) an index for the size class at the end of the year; and (4) the number
    of individuals that made that transition.  Deaths are specified by a zero in 
    the third column.}
\item{\code{fec}}{A matrix specifing the fertility data.  It has four columns:
    (1) the year; (2) an index for the size class of the parent; (3) the number
    of parents in that class; and (4) the number of newborns coming from the 
    parents in that class.}
\item{\code{p0}}{Survival of newborns.  Must be specified if \code{census = "pre"}.}
\item{\code{classes}}{A character vector with the name of each age class.  Must be 
    specified.}
\item{\code{census}}{Whether to construct a model that represents a pre-breeding census
    or a post-breeding census.  Defaults to \code{"pre"}.}
\item{\code{type}}{Type of matrix.  Must be set to \code{"size"}.}
\item{\code{data.type}}{Must be set to \code{"aggregated"}.}
\item{\code{fec.sampling}}{Whether the distribution of the number of offspring among
    mothers should be assumed to be Poisson or binomial.  Use the latter if at
    most one offspring is produced in a year.  Defaults to \code{"Poisson"}.}
\item{\code{seedbank}}{Logical.  Indicates whether there is a persistent seed bank.  If 
    \code{TRUE}, then values must be specified for \code{seed.germ} and 
    \code{seed.surv}.  Defaults to \code{FALSE}.}
\item{\code{seed.germ}}{The fraction of surviving seeds in the seedbank that germinate.
    May be omitted if \code{seedbank = FALSE}.}
\item{\code{seed.surv}}{The annual survival of seeds in the seedbank.
    May be omitted if \code{seedbank = FALSE}.}
\item{\code{annual}}{Logical.  Whether to extract annual models as well as the mean model.
    If set to \code{TRUE}, then \code{years} must be specified.}
\item{\code{years}}{A numeric list of the years for which annual models should be built.
    These must match values in the first column of \code{x}.}
}

This is somewhat complicated.  See the "Examples" for insights on how to actually
do it.  Examples of the data structure for the aggregated option are in the datasets
\code{Cereus.Grow} and \code{Cereus.Fruit}.
}

\value{
An object of class dmm.  This is a list containing the following components:

\code{mat}        The transition matrix

\code{classes}    The class names

\code{type}       The type of model (age or size)

\code{census}     The type of census (pre or post)

\code{surv}       Survival of each class

\code{fec}        Fertility of each class

\code{grow}       Growth matrix

\code{grow.binom} Growth parameters, set up as nested binomials

\code{se}         Standard errors of the parameter estimates, itself a list with
    components surv, fec, and grow (the latter are the SE of the grow.binom
    parameters)
    
\code{mat.ann}    A list of the annual matrices, each of which is itself a dmm
    object

Some comments:
\itemize{
\item{The classes, surv, and fec components include newborns, even if the model
    is for a pre-breeding census}
\item{The grow and grow.binom components are NULL for age-structured models}
\item{The se component is non-NULL only when data.type="aggregated"}
}
}


\note{
This function is something of a kludge, with different pieces added hastily
on an ad hoc basis.  A lot of potentially desireable combinations of options are
unavailable.  This will be rationalized in a future version.
}

\author{Bruce E. Kendall (\email{kendall@bren.ucsb.edu})}

\seealso{
A print method is specified in \code{\link{print.dmm}}. \code{\link{dmm.iter}} and 
\code{\link{dmm.asymp}} perform deterministic analyses of demographic
matrix models.  \code{\link{dmm.PVA}} and \code{\link{tulja.dmm}} perform
stochastic analyses of matrix models.
}

\examples{
## Construct an age-structured model
data(sealion)
sealion.dmm <- dmm(sealion, type="age", census="pre", terminal=F)

## Construct a size structured model from rate data.  
## Read in the Hudsonia data
data(Hudsonia.Grow)
data(Hudsonia.Seeds)

# Choose some class boundaries
classmax <- c(10,20,30,40,60,100,150,500)
classmin <- c(0,10,20,30,40,60,100,150)

## Convert the sizes into classes
Hudsonia.g2 <- dmm.sizeclass(Hudsonia.Grow,classmin,classmax)

## Look at size-dependent survival
Hudsonia.surv <- dmm.survplot(Hudsonia.Grow,classmin,classmax)

## Look at the plot of seeds vs. area, and do the regression
plot(Hudsonia.Seeds$Area,sqrt(Hudsonia.Seeds$Seeds), xlab="Area", ylab="sqrt(Number of seeds)")
fec.lm <- lm(sqrt(Seeds) ~ Area, data=Hudsonia.Seeds )

## Calculate the seed production for the class medians
fec.med <- predict(fec.lm,newdata=data.frame(Area=Hudsonia.surv$median.size))^2

# Build the lifetable
Hudsonia.lt <- cbind(2:9,Hudsonia.surv$survival,fec.med)

# Add the seedlings
Hudsonia.lt <- rbind(c(1,(0.5+0.4545)/2,0),Hudsonia.lt)

# Add the seedlings to the growth matrix
Hudsonia.gmat <- Hudsonia.g2$gmat
Hudsonia.gmat <- rbind(rep(0,8),Hudsonia.gmat)
Hudsonia.gmat <- cbind(c(0,1,rep(0,7)),Hudsonia.gmat)

# Finally, construct the model
Hudsonia.dmm <- dmm(Hudsonia.lt, grow=Hudsonia.gmat, type="size",
    census="pre", data.type="rates", seedbank=T, seed.germ=0.01, seed.surv=0.5,
    classes=c("seedBnk","seedlng","teensy","tiny","little","small","medium",
                "large","huge","giant"))


## Read in the Cereus data
data(Cereus.Grow)
data(Cereus.Fruit)

## Construct size structured model using aggregated data, and get annual models
cereus2.dmm <- dmm(Cereus.Grow, fec=Cereus.Fruit, p0=0.03448, type="size",
    census="pre", data.type="aggregated",
    classes=c("<1","1-4","4-8","8-16","16-32",">32"),annual=T,years=1988:1992)

## Look at the annual matrices
cereus2.dmm$mat.ann

## Extract the first annual matrix as an independent dmm object
A1.dmm <- cereus2.dmm$mat.ann[[1]]
}

\keyword{classes}

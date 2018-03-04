#' Estimate SEG model parameters
#' 
#' Given a time series of abundances, this function gives estimates and 
#' confidence intervals for the parameters of the stochastic exponential growth 
#' (SEG) parameters. It uses the "linear regression approach" described on pp. 
#' 120-121 of Dennis et al. (1991) (see also Chapter 3 of Morris & Doak 2002 for
#' a friendlier description). By default, the time series of abundances is 
#' assumed to be equally spaced in time, with a timestep of one time unit; 
#' deviations from this can be accomodated by providing a vector of times 
#' associated with each abundance estimate.
#' 
#' 
#' @param Nt Time series of abundance. Must be a numeric vector of positive
#'   values, with length at least 3
#' @param Yeart Sequence of times ("years") at which abundance is estimated
#'   (must be the same length as `Nt`). If the abundance estimates are available
#'   at each time step, this can be omitted
#' @param alpha Significance level for confidence levels (default is 0.05).
#'   
#' @return A list with the following elements:
#'   
#'   \item{mu}{Estimate of the mean growth rate}
#'   
#'   \item{sigma2}{Estimate of the variance in the growth rate}
#'   
#'   \item{params.CI}{A matrix with the CI of mu in the first row and the CI of 
#'   sigma2 in the second row}
#'   
#'   \item{alpha}{The alpha level passed in by the user}
#' 
#' @references Dennis, B., Munholland, P. L. and Scott, J. M. (1991) Estimation 
#'   of growth and extinction parameters for endangered species. 
#'   \emph{Ecological  Monographs} \bold{61}, 115--143.
#'   
#'   Morris, W. F. and Doak D. F. (2002) \emph{Quantitative Conservation 
#'   Biology: Theory and Practice of Population Viability Analysis}. Sunderland:
#'   Sinauer Associates.
#'   
#' @examples
#' data(grizzly)
#' estimate_SEG_params(grizzly$N)
#'   
estimate_SEG_params <- function(Nt, Yeart = seq(along = Nt), alpha = 0.05) {
  # Error checking
  stopifnot(length(Nt) == length(Yeart),
            alpha > 0,
            alpha < 1)
  
  # Calculate the elements needed for the regression
  r <- diff(log(Nt))
  x <- sqrt(diff(Yeart))
  y <- r / x
  
  # Do the regression and extract the components
  xy.lm <- lm(y ~ x - 1) # The '-1' specifies no intercept
  mu <- xy.lm$coef[1]
  sigma2 <- var(xy.lm$resid)
  param_CI <- confint(xy.lm, "x", 1 - alpha)
  DF <- length(r) - 1
  param_CI <- rbind(param_CI,
                    DF * sigma2 / qchisq(c(1 - alpha/2, alpha/2), DF))
  rownames(param_CI) <- c("mu", "sigma2")
  return(list(mu = mu, sigma2 = sigma2, params.CI = param_CI, alpha = alpha))
}
# Define the function that calculates the CDF of extinction risk
# There is a pinvGauss in the SuppDists package, but it has a different
#   parameterization and only works for negative values of mu

pinvGauss <- function(T,mu,s2,d)
    pnorm((-d - mu*T)/sqrt(s2*T)) + exp(-2*mu*d/s2) * pnorm((-d + mu*T)/sqrt(s2*T))

AICc <- function(x)
{
    if (class(x) == "lm") {
        p <- length(coef(x)) + 1
        Q <- length(fitted(x))
        AIC(x) - 2*p + 2*p*Q/(Q-p-1)
    } else if (class(x) == "nls") {
        p <- length(coef(x)) + 1
        Q <- length(fitted(x))
        AIC(x) - 2*p + 2*p*Q/(Q-p-1) + 2
    } else {
        cat("AICc is not implemented for objects of class",class(x),"\n")
#        NULL
    }
}

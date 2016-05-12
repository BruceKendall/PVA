bc <- function(x,theta) {
    if (theta == 0) {
        log(x)
    } else {
        x^theta
    }
}

profile.theta <- function(rt,Nt) {
    i <- 1
    fit <- vector("numeric",51)
    for (theta in seq(-2.5,2.5,.1)) {
        N <- bc(Nt,theta)
        test.lm <- lm(rt~N)
        fit[i] <- var(resid(test.lm))
        i <- i+1
    }
    plot(seq(-2.5,2.5,.1),fit,type='b',xlab="theta",ylab="residual variance")
}

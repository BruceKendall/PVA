dmm.iter <- function(x, N0, steps, col=1:6, ...)
{
    if (class(x) == "dmm") {
        mat <- x$mat
        classes <- x$classes
    } else if (class(x) == "matrix") {
        mat <- x
        classes <- 1:(dim(mat)[1])
    } else {
        stop("x must be of class dmm or matrix")
    }
    N <- N0
    N <- rbind(N,t(mat%*%N))
    for (i in 2:steps) 
        N <- rbind(N,t(mat%*%N[i,]))
    
    Year <- 0:steps
    
    matplot(Year,N,type='l',col=col,...)
    legend(steps*.75,min(N)+diff(range(N))/2,as.character(classes),col=col,
        lty=1:min(c(5,length(N0))))
    invisible(N)
}

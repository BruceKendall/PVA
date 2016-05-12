plot.ext.risk <- function(x, ...) {
    if (is.null(x$alpha)) {
        plot(x$year, x$extrisk, type='l', xlab="Years", ylab="Cumulative extinction risk",...)
    } else {
        last <- length(x$year)
        par(mfrow=c(2,2))
        matplot(x$ext.ci, type='l', lty=2, xlab="Years", ylab="Cumulative extinction risk",col=1,
                main=paste("MLE and ",(1-x$alpha)*100,"% CI"), ...)
        lines(x$extrisk)
        matplot(t(apply(x$ext.boot,1,quantile,c(0.25,.5, 0.55),na.rm=T)),type='l', lty=1, xlab="Years", 
                ylab="Cumulative extinction risk",col=1, main="Median and quartiles",...)
        use <- !is.na(x$ext.boot[last,])
        dens <- density(x$ext.boot[last,use])
        plot(dens,main="Extinction risk at end time", xlab="Extinction risk",
                ylab="PDF")
        plot(sort(x$ext.boot[last,use]),(1:(sum(use)))/sum(use),main="Extinction risk at end time", 
                type='l', xlab="Extinction risk", ylab="CDF")
        par(mfrow=c(1,1))
    }
}

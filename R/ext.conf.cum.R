ext.conf.cum <- function(x, year=length(x$year)) 
{
    use <- !is.na(x$ext.boot[year,])
    data.frame(Extinction.risk=c(0,sort(x$ext.boot[year,use])),
        Cumulative.confidence=c(0,(1:(sum(use)))/sum(use)))
}

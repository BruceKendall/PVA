dmm.growthplot <- function(x, ...)
{
    nyr <- dim(x)[2] - 1
    grow <- cbind(x[,2],x[,3])
    if (nyr > 2)
        for (i in 3:nyr)
            grow <-rbind(grow, cbind(x[,i],x[,i+1]))
    # get rid of deaders
    grow <- grow[grow[,1]*grow[,2] > 0,]
    
    plot(grow[,1],grow[,2],tck=1,xlab="Size in year t", ylab="Size in year t+1", ...)
    invisible(grow)
}

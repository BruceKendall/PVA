dmm.survplot <- function(x, classmin=NULL, classmax=NULL, ...)
{
    nyr <- dim(x)[2] - 1
    grow <- cbind(x[,2],x[,3])
    if (nyr > 2)
        for (i in 3:nyr)
            grow <-rbind(grow, cbind(x[,i],x[,i+1]))
    # get rid of deaders
    grow <- grow[grow[,1] > 0,]
    # Now recode for survival
    grow[grow[,2] > 0,2] <- 1
    
    plot(grow[,1],grow[,2], xlab="Size", ylab="Survival",pch='|',cex=.5,...)
    grow <- data.frame(size=grow[,1],surv=grow[,2])
    survmod <- glm(surv~size, data=grow, family=binomial)
    vec <- order(grow[,1])
    lines(grow[vec,1],survmod$fitted[vec])
    survmod.lo <- loess(grow[vec,2]~grow[vec,1],span=.75)
    lines(grow[vec,1],survmod.lo$fitted,lty=2)

    med <- vector("numeric",length(classmin))    
    smed <- vector("numeric",length(classmin))    
    if (!is.null(classmin)) {
        for (i in 1:length(classmin)) {
            vec <- grow[,1] > classmin[i] & grow[,1] < classmax[i]
            med[i] <- median(grow[vec,1])
            smed[i] <- predict(survmod,newdata=data.frame(size=med[i]),type="response")
            lines(c(classmin[i],classmax[i]),c(smed[i],smed[i]),col='red',lwd=2)
            if (i > 1) lines(c(classmin[i],classmin[i]),c(smed[i-1],smed[i]),col='red',lwd=2)
        }
        data.frame(median.size=med,survival=smed)
    }
}

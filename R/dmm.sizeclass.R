dmm.sizeclass <- function(x,classmin,classmax)
{
    nyr <- dim(x)[2] - 1
    nind <- dim(x)[1]
    nclass <- length(classmin)
    grow <- cbind(rep(1,nind),x[,2],x[,3])
    if (nyr > 2)
        for (i in 3:nyr)
            grow <-rbind(grow, cbind(rep(i-1,nind),x[,i],x[,i+1]))
    # get rid of deaders
    grow <- grow[grow[,2] > 0,]
    nyr <- nyr-1
    growc <- matrix(0,nclass*(nclass+1)*nyr,4)
    for (i in 1:nyr) {
        growi <- grow[grow[,1]==i,]
        for (j in 1:nclass) {
            vecj <- growi[,2]>classmin[j] & growi[,2]<classmax[j]
            if (sum(vecj) > 0) {
                growj <- growi[vecj,]
                if (is.vector(growj)) growj <- t(as.matrix(growj))
                growc[(i-1)*nclass*(nclass+1)+(j-1)*(nclass+1)+1,] <- c(i, j, 0, 
                        sum(growj[,3]==0))
                for (k in 1:nclass) {
                    growc[(i-1)*nclass*(nclass+1)+(j-1)*(nclass+1)+k+1,] <- c(i, j, k, 
                        sum(growj[,3]>classmin[k] & growj[,3]<classmax[k]))
                }
            }
        }
    }
    x <- growc
    from.class <- sort(unique(x[,2]))
    nclass <- length(from.class)
    gbar <- matrix(0,nclass,nclass)
    N <- vector("numeric",nclass)
    for (i in from.class) {
        ivec <- x[,2] == i
        jvec <- x[,3] > 0
        N[i] <- sum(x[ivec,4])
        Ns <- sum(x[ivec&jvec,4])
        for (j in from.class) {
            ijvec <- ivec & x[,3]==j
            Ng <- sum(x[ijvec,4])
            gbar[j,i] <- Ng/Ns
        }
    }

    list(growc=growc, gmat=gbar)
}

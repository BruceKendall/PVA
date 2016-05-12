VC.dmm <- function(x)
{
    nmat <- length(x$mat.ann)
    nfec <- length(x$fec)
    nsurv <- length(x$surv)
    nclass <- length(x$classes)
    ngrow <- nclass*(nclass-1)
    
    fecmat <- NULL
    survmat <- NULL
    growmat <- NULL
    for (i in 1:nmat) {
        A <- x$mat.ann[[i]]
        fecmat <- cbind(fecmat,A$fec)
        survmat <- cbind(survmat,A$surv)
        growmat <- cbind(growmat,as.vector(A$grow.binom))
    }
    parmat <- rbind(fecmat,survmat,growmat)
    pbar <- apply(parmat,1,mean)
    npar <- dim(parmat)[1]
    VC <- matrix(0,npar,npar)
    cormat <- matrix(0, npar, npar)
    parmat[is.nan(parmat)] <- 0
    options(warn=-1)
    for (i in 1:npar) {
        p1 <- parmat[i,]
        for (j in 1:npar) {
            p2 <- parmat[j,]
            VC[i,j] <- cov(p1,p2)
            cormat[i,j] <- cor(p1,p2)
        }
    }
    options(warn=0)
    cormat[is.na(cormat)] <- 0
    list(VC=VC, cormat=cormat, pbar=pbar)
}

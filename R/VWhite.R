VWhite <- function(x) {
    nmat <- length(x$mat.ann)
    nfec <- length(x$fec)
    nsurv <- length(x$surv)
    nclass <- length(x$classes)
    ngrow <- nclass^2
    
    fecmat <- NULL
    survmat <- NULL
    growmat <- NULL
    fecse <- NULL
    survse <- NULL
    growse <- NULL
    for (i in 1:nmat) {
        A <- x$mat.ann[[i]]
        fecmat <- cbind(fecmat,A$fec)
        survmat <- cbind(survmat,A$surv)
        growmat <- cbind(growmat,as.vector(A$grow.binom))
        fecse <- cbind(fecse,A$se$fec)
        survse <- cbind(survse,A$se$surv)
        growse <- cbind(growse,as.vector(A$se$grow))
    }
    growmat[is.nan(growmat)] <- 0
    growse[is.na(growse)] <- 0
    VCfec <- vector("numeric",nfec)
    Vfec <- apply(fecmat,1,var)
    VCsurv <- vector("numeric",nsurv)
    Vsurv <- apply(survmat,1,var)
    VCgrow <- vector("numeric",ngrow)
    Vgrow <- apply(growmat,1,var)
    
    for (i in 1:nfec) {
        fbar <- mean(fecmat[i,])
        if (Vfec[i] > 0) {
            VC <- optimize(function(x,fbar,fi,Vi){(sum((fi-fbar)^2/(x+Vi)) - (length(fi)-1))^2}, 
                c(0,2*Vfec[i]),
                fbar=fbar,fi=fecmat[i,],Vi=fecse[i,]^2)
            VCfec[i] <- VC$min
        } else {
            VCfec[i] <- 0
        }
    }
    
    for (i in 1:nsurv) {
        fbar <- mean(survmat[i,])
        if (Vsurv[i] > 0) {
            VC <- optimize(function(x,fbar,fi,Vi){(sum((fi-fbar)^2/(x+Vi)) - (length(fi)-1))^2}, 
                c(0,2*Vsurv[i]),
                fbar=fbar,fi=survmat[i,],Vi=survse[i,]^2)
            VCsurv[i] <- VC$min
        } else {
            VCsurv[i] <- 0
        }
    }
    
    for (i in 1:ngrow) {
        fbar <- mean(growmat[i,])
        if (Vgrow[i] > 0) {
            VC <- optimize(function(x,fbar,fi,Vi){(sum((fi-fbar)^2/(x+Vi)) - (length(fi)-1))^2}, 
                c(0,2*Vgrow[i]),
                fbar=fbar,fi=growmat[i,],Vi=growse[i,]^2)
            VCgrow[i] <- VC$min
        } else {
            VCgrow[i] <- 0
        }
    }
    list(VCfec=VCfec,Vfec=Vfec,VCsurv=VCsurv,Vsurv=Vsurv,VCgrow=VCgrow,Vgrow=Vgrow)
}

sens.vr <- function(x) {
    sens.asymp <- dmm.asymp(x, quiet=T)
    S <- sens.asymp$sensmat
    lam1 <- sens.asymp$lambda1
    nclass <- length(x$classes)
    if (x$type == "age") {
        if (x$census == "pre" ) {
            S.fec <- vector("numeric",nclass+1)
            S.surv <- vector("numeric",nclass+1)
            S.fec[1] <- 0
            S.surv[1] <- sum(x$fec[-1] * S[1,])
            for (i in 1:(nclass-1)) {
                S.fec[i+1] <- x$surv[1] * S[1,i]
                S.surv[i+1] <- S[i+1,i]
            }
            S.fec[nclass+1] <- x$surv[1] * S[1,nclass]
            S.surv[nclass+1] <- 0
            if (x$mat[nclass,nclass] > 0) S.surv[nclass+1] <- S[nclass,nclass]
        } else {
            S.fec <- vector("numeric",nclass)
            S.surv <- vector("numeric",nclass)
            S.fec[1] <- 0
            S.surv[1] <- S[2,1] + x$fec[2] * S[1,1]
            for (i in 2:(nclass-1)) {
                S.fec[i] <- x$surv[i-1] * S[1,i-1]
                S.surv[i] <- S[i+1,i] + x$fec[i+1] * S[1,i]
            }
            if (x$mat[nclass,nclass] > 0) {
                S.surv[nclass] <- S[nclass,nclass] + x$fec[nclass] * S[1,nclass]
                S.fec[nclass] <- x$surv[nclass-1] * S[1,nclass-1] + 
                    x$surv[nclass] * S[1,nclass]
            } else {
                S.surv[nclass] <- 0
                S.fec[nclass] <- x$surv[nclass-1] * S[1,nclass-1] 
            }
        }
        E.fec <- S.fec * x$fec / lam1
        E.surv <- S.surv * x$surv / lam1
        S.grow <- NULL
        E.grow <- NULL
        S.grow.binom <- NULL
        E.grow.binom <- NULL
    } else {
        if (x$census == "pre" ) {
            S.fec <- vector("numeric",nclass+1)
            S.surv <- vector("numeric",nclass+1)
            S.grow <- matrix(0, nclass, nclass)
            S.grow.binom <- matrix(0, nclass, nclass)
            S.fec[1] <- 0
            S.surv[1] <- sum(x$fec[-1] * S[1,])
            for (i in 1:nclass) {
                S.fec[i+1] <- x$surv[1] * S[1,i]
                S.surv[i+1] <- sum(S[,i] * x$grow[,i])
                S.grow[,i] <- x$surv[i+1] * S[,i]
                vec <- x$grow.binom[,i] > 0 & x$grow.binom[,i] < (1-.00001)
                jvec <- (1:nclass)[vec]
                jvec <- min(jvec):max(jvec)
                gvec <- x$grow.binom[jvec,i]
                ng <- length(gvec)
                for (j in 1:ng) {
                    if (x$grow.binom[jvec[j],i] > 0) {
                        S.grow.binom[jvec[j],i] <- S[jvec[j],i] * x$surv[i+1] *
                            x$grow[jvec[j],i]/x$grow.binom[jvec[j],i]
                    } else {
                        S.grow.binom[jvec[j],i] <- 0
                    }
                    for (k in j:ng) {
                        S.grow.binom[jvec[j],i] <- S.grow.binom[jvec[j],i] -
                            S[jvec[k]+1,i] * x$surv[i+1] *
                            x$grow[jvec[k]+1,i]/(1 - x$grow.binom[jvec[j],i])
                    }
                }
            }
        } else {
            S.fec <- vector("numeric",nclass)
            S.surv <- vector("numeric",nclass)
            S.grow <- matrix(0, nclass, nclass)
            S.grow.binom <- matrix(0, nclass, nclass)
            for (i in 1:nclass) {
                S.fec[i] <- sum(x$surv * x$grow[i,] * S[1,i])
                S.surv[i] <- S[,i]* x$grow[,1] + sum(S[1,] * x$fec * x$grow[,i])
                S.grow[,i] <- x$surv[i] * S[,i] + x$surv[i] * fec * S[1,]
                vec <- x$grow.binom[,i] > 0 & x$grow.binom[,i] < (1-.00001)
                jvec <- (1:nclass)[vec]
                jvec <- min(jvec):max(jvec)
                gvec <- x$grow.binom[jvec,i]
                ng <- length(gvec)
                for (j in 1:ng) {
                    if (x$grow.binom[jvec[j],i] > 0) {
                        S.grow.binom[jvec[j],i] <- 
                            (S[jvec[j],i] + fec[jvec[j]] * S[1,jvec[j]]) * x$surv[i+1] *
                            x$grow[jvec[j],i]/x$grow.binom[jvec[j],i]
                    } else {
                        S.grow.binom[jvec[j],i] <- 0
                    }
                    for (k in j:ng) {
                        S.grow.binom[jvec[j],i] <- S.grow.binom[jvec[j],i] -
                            (S[jvec[k]+1,i] + fec[jvec[k]+1] * S[1,jvec[k]+1]) * x$surv[i+1] *
                            x$grow[jvec[k]+1,i]/(1 - x$grow.binom[jvec[j],i])
                    }
                }
            }
        }
        E.fec <- S.fec * x$fec / lam1
        E.surv <- S.surv * x$surv / lam1
        E.grow <- S.grow * x$grow / lam1
        E.grow.binom <- S.grow.binom * x$grow.binom / lam1
    }
    list(S.fec=S.fec,S.surv=S.surv,S.grow=S.grow,S.grow.binom=S.grow.binom,
        E.fec=E.fec,E.surv=E.surv,E.grow=E.grow,E.grow.binom=E.grow.binom)
}

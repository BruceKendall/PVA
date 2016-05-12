tulja.dmm <- function(x)
{
    lam1 <- dmm.asymp(x,quiet=T)$lambda1
    VC <- VC.dmm(x)$VC
    S.vr <- sens.vr(x)
    S <- c(S.vr$S.fec,S.vr$S.surv,as.vector(S.vr$S.grow.binom))
    tau2 <- 0
    npar <- length(S)
    for (i in 1:npar) 
        for (j in 1:npar) 
            tau2 <- tau2 + S[i]*S[j]*VC[i,j]
    log(lam1) - 0.5 * tau2/lam1^2
}

count.DI.PVA <- function(Nt, Yeart=seq(along=Nt), Nc=Nt[length(Nt)], Nx=1, 
    year.max=100, nboot=0, alpha=0.05, meas.err=c("none","simple"), Nt.se=NULL, ...) {
    n = length(Nt)
    m <- n-1         # m is q from Morris and Doak
    d <- log(Nc) - log(Nx)
    rt = log(Nt[-1]/Nt[-n])
    if (meas.err[1] == "simple")
        rt.se <- (Nt.se[-1]/Nt[-1])^2 + (Nt.se[-n]/Nt[-n])^2 
    if (n == diff(range(Yeart))+1) {
        # Calculate mu and sigma2 directly
        mu <- mean(rt)
        sigma2 <- var(rt)
        if (meas.err[1] == "simple") 
            sigma2 <- max(c(0.001, sigma2 - mean(rt.se)))
    } else {
        ### Estimation using regression
        # Create x and y
        x <- sqrt(Yeart[2:(m+1)]-Yeart[1:m])
        y <- rt/x
        # Do the regression and extract the components
        xy.lm <- lm(y~x-1) # The '-1' specifies no intercept
        mu <- xy.lm$coef[1]
        sigma2 <- var(xy.lm$resid)
    }

    # Calculate the CDF of extinction
    extrisk <- pinvGauss(1:year.max,mu,sigma2,d)
    plot(1:year.max, extrisk, type='l', xlab="Years", ylab="Cumulative extinction risk",...)

    if (nboot > 0) {
        if (n == diff(range(Yeart))+1) {
            iboot = matrix(sample(m,nboot*m,T),m,nboot)
            rboot = matrix(rt[iboot],m,nboot)
            mu.boot = apply(rboot, 2, mean)
            s2.boot = apply(rboot, 2, var)
            if (meas.err[1] == "simple") {
                r.se.boot = matrix(rt.se[iboot],m,nboot)
                s2.boot <- s2.boot - apply(r.se.boot, 2, mean)
                s2.boot[s2.boot <= 0.001] <- 0.001
            }
        } else {
            iboot = matrix(sample(m,nboot*m,T),m,nboot)
            x.boot = matrix(x[iboot],m,nboot)
            y.boot = matrix(y[iboot],m,nboot)
            mu.boot <- vector("numeric",nboot)
            s2.boot <- vector("numeric",nboot)
            for (k in 1:nboot) {
                boot.lm <- lm(y.boot[,k]~x.boot[,k]-1)
                mu.boot[k] <- boot.lm$coef[1]
                s2.boot[k] <- var(boot.lm$resid)
            }
        }
        TT <- matrix(1:year.max,year.max,nboot)         # Each column of these matrices
        mm <- matrix(mu.boot,year.max,nboot,byrow=T)    #   is a separate replicate; the
        ss <- matrix(s2.boot,year.max,nboot,byrow=T)    #   rows are time into future
        print(c(sum(is.na(mm)),sum(is.na(ss))))
        print(c(range(mm),range(ss)))
        ext.boot <- pinvGauss(TT,mm,ss,d)
        print(sum(is.na(ext.boot)))

        # Calculate and plot the CI
        ext.ci <- apply(ext.boot,1,quantile,c(alpha/2,1-alpha/2),na.rm=TRUE)
        matplot(t(ext.ci), type='l', lty=2, xlab="Years", ylab="Cumulative extinction risk",ylim=c(0.00000001,max(ext.boot,na.rm=T)),col=1,...)
        lines(extrisk)
        output <- list(extrisk=extrisk,year=1:year.max,nboot=nboot,ext.boot=ext.boot,alpha=alpha,ext.ci=t(ext.ci),
                        model.type="Density independent count",model.params=c(mu=mu,sigma2=sigma2))
        class(output) <- "ext.risk"
        output
    } else {
        output <- list(extrisk=extrisk,year=1:year.max,nboot=NULL,ext.boot=NULL,alpha=NULL,ext.ci=NULL,
                        model.type="Density independent count",model.params=c(mu=mu,sigma2=sigma2))
        class(output) <- "ext.risk"
        output
    }
}

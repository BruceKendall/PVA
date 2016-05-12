count.DD.PVA <- function(Nt, Yeart=seq(along=Nt), Nc=Nt[length(Nt)], Nx=1, 
    year.max=100, nsim=1000, nboot=0, alpha=0.05, 
    boot.type=c("nonparametric","parameteric"),
    model=c("Ricker","theta","logistic","Beverton","delay","DI"), 
    enoise=c("normal","resample","DS"),
    par.start=NULL, meas.err=FALSE, ...) 
{
    model <- model[1]
    enoise <- enoise[1]
    boot.type <- boot.type[1]
    Q = length(Nt)-1
    rt <- diff(log(Nt))
    if (model == "delay" && length(Nc) == 1) Nc <- c(Nt[Q],Nt[Q+1])
    Nt <- Nt[1:Q]
    if (boot.type == "parametric") require(MASS)

    if (model == "DI") {
        fit.model <- lm(rt ~ 1)
        mu <- coef(fit.model)[1]
        s2 <- var(resid(fit.model))
        if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
        sims <- sim.DI(Nc, c(mu, s2), year.max, nsim, enoise, resid(fit.model))
    } else if (model == "Ricker") {
        fit.model <- lm(rt ~ Nt)
        mu <- coef(fit.model)[1]
        K <- -mu / coef(fit.model)[2]
        s2 <- var(resid(fit.model))
        if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
        sims <- sim.Ricker(Nc, c(mu, K, s2), year.max, nsim, enoise, resid(fit.model))
    } else if (model == "theta") {
        if (is.null(par.start)) { #Try starting from the Ricker values
            fit.model <- lm(rt ~ Nt)
            mu <- coef(fit.model)[1]
            K <- -mu / coef(fit.model)[2]
            par.start <- list(mu=mu, K=K, theta=1)
        } 
        fit.model <- nls(rt ~ mu*(1-(Nt/K)^theta), start=par.start, ...)
        mu <- coef(fit.model)[1]
        K <- coef(fit.model)[2]
        theta <- coef(fit.model)[3]
        s2 <- var(resid(fit.model))
        if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
        sims <- sim.theta(Nc, c(mu, K, theta, s2), year.max, nsim, enoise, resid(fit.model))
    } else if (model == "logistic") {
        if (is.null(par.start)) { #Try starting from the Ricker values
            fit.model <- lm(rt ~ Nt)
            mu <- coef(fit.model)[1]
            K <- -mu / coef(fit.model)[2]
            par.start <- list(R=exp(mu), K=K)
        } 
        fit.model <- nls(rt ~ log(R*(1-(Nt/K))), start=par.start, ...)
        R <- coef(fit.model)[1]
        K <- coef(fit.model)[2]
        s2 <- var(resid(fit.model))
        if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
        sims <- sim.logistic(Nc, c(R, K, s2), year.max, nsim, enoise, resid(fit.model))
    } else if (model == "Beverton") {
        if (is.null(par.start)) { #Try starting from the Ricker values
            fit.model <- lm(rt ~ Nt)
            mu <- coef(fit.model)[1]
            K <- -mu / coef(fit.model)[2]
            par.start <- list(R=exp(mu), a=(exp(mu)-1)/K)
        } 
        fit.model <- nls(rt ~ log(R/(1+a*Nt)), start=par.start, ...)
        R <- coef(fit.model)[1]
        a <- coef(fit.model)[2]
        s2 <- var(resid(fit.model))
        if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
        sims <- sim.Beverton(Nc, c(R, a, s2), year.max, nsim, enoise, resid(fit.model))
    } else if (model == "delay") {
        fit.model <- lm(rt[2:Q] ~ Nt[2:Q] + Nt[1:(Q-1)])
        s2 <- var(resid(fit.model))
        if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
        sims <- sim.delay(Nc, c(coef(fit.model), s2), year.max, nsim, enoise, resid(fit.model))
    } else {
        stop("Model of type '",model,"' is not supported in count.DD.PVA.", call.=F)
    }

    # Calculate the CDF of extinction
    sims[is.nan(sims)] <- 0
    sims[is.infinite(sims)] <- 0
    minsize <- apply(sims,2,cummin)
    QE <- (minsize < Nx)
    extrisk <- apply(QE,1,sum)/nsim
    plot(1:year.max, extrisk, type='l', xlab="Years", ylab="Cumulative extinction risk",...)

    parhat <- coef(fit.model)
    if (model=="Ricker") { 
        parhat <- c(mu,K)
        names(parhat) <- c("mu","K")
    }
    
    if (nboot > 0 && boot.type =="nonparametric") {
        iboot = matrix(sample(Q,nboot*Q,T),Q,nboot)
        Nt.boot = matrix(Nt[iboot],Q,nboot)
        rt.boot = matrix(rt[iboot],Q,nboot)
        ext.boot <- matrix(0,year.max,nboot)
        if (model == "DI") {
            for (i in 1:nboot) {
                fit.model <- lm(rt.boot[,i] ~ 1)
                mu <- coef(fit.model)[1]
                s2 <- var(resid(fit.model))
                if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
                sims <- sim.DI(Nc, c(mu, s2), year.max, nsim, enoise, resid(fit.model))
                sims[is.nan(sims)] <- 0
                sims[is.infinite(sims)] <- 0
                minsize <- apply(sims,2,cummin)
                QE <- (minsize < Nx)
                ext.boot[,i] <- apply(QE,1,sum)/nsim
            }
        } else if (model == "Ricker") {
            for (i in 1:nboot) {
                fit.model <- lm(rt.boot[,i] ~ Nt.boot[,i])
                mu <- coef(fit.model)[1]
                K <- -mu / coef(fit.model)[2]
                s2 <- var(resid(fit.model))
                if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
                sims <- sim.Ricker(Nc, c(mu, K, s2), year.max, nsim, enoise, resid(fit.model))
                sims[is.nan(sims)] <- 0
                sims[is.infinite(sims)] <- 0
                minsize <- apply(sims,2,cummin)
                QE <- (minsize < Nx)
                ext.boot[,i] <- apply(QE,1,sum)/nsim
            }
        } else if (model == "theta") {
            for (i in 1:nboot) {
                fit.model <- nls(rt.boot[,i] ~ mu*(1-(Nt.boot[,i]/K)^theta), start=parhat, ...)
                mu <- coef(fit.model)[1]
                K <- coef(fit.model)[2]
                theta <- coef(fit.model)[3]
                s2 <- var(resid(fit.model))
                if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
                sims <- sim.theta(Nc, c(mu, K, theta, s2), year.max, nsim, enoise, resid(fit.model))
                sims[is.nan(sims)] <- 0
                sims[is.infinite(sims)] <- 0
                minsize <- apply(sims,2,cummin)
                QE <- (minsize < Nx)
                ext.boot[,i] <- apply(QE,1,sum)/nsim
            }
        } else if (model == "logistic") {
            for (i in 1:nboot) {
                fit.model <- nls(rt.boot[,i] ~ log(R*(1-(Nt.boot[,i]/K))), start=parhat, ...)
                R <- coef(fit.model)[1]
                K <- coef(fit.model)[2]
                s2 <- var(resid(fit.model))
                if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
                sims <- sim.logistic(Nc, c(R, K, s2), year.max, nsim, enoise, resid(fit.model))
                sims[is.nan(sims)] <- 0
                sims[is.infinite(sims)] <- 0
                minsize <- apply(sims,2,cummin)
                QE <- (minsize < Nx)
                ext.boot[,i] <- apply(QE,1,sum)/nsim
            }
        } else if (model == "Beverton") {
            for (i in 1:nboot) {
                fit.model <- nls(rt.boot[,i] ~ log(R/(1+a*Nt.boot[,i])), start=parhat, ...)
                R <- coef(fit.model)[1]
                a <- coef(fit.model)[2]
                s2 <- var(resid(fit.model))
                if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
                sims <- sim.Beverton(Nc, c(R, a, s2), year.max, nsim, enoise, resid(fit.model))
                sims[is.nan(sims)] <- 0
                sims[is.infinite(sims)] <- 0
                minsize <- apply(sims,2,cummin)
                QE <- (minsize < Nx)
                ext.boot[,i] <- apply(QE,1,sum)/nsim
            }
        } else if (model == "delay") {
            iboot = matrix(sample(2:Q,nboot*(Q-1),T),Q-1,nboot)
            Nt.boot = matrix(Nt[iboot],Q-1,nboot)
            Ntm1.boot = matrix(Nt[iboot-1],Q-1,nboot)
            rt.boot = matrix(rt[iboot],Q-1,nboot)
            ext.boot <- matrix(0,year.max,nboot)
            for (i in 1:nboot) {
                fit.model <- lm(rt.boot[,i] ~ Nt.boot[,i]+Ntm1.boot[,i])
                s2 <- var(resid(fit.model))
                if (enoise == "DS") s2 <- as.vector(coef(lm(resid(fit.model)^2~I(1/Nt))))
                sims <- sim.delay(Nc, c(coef(fit.model), s2), year.max, nsim, enoise, resid(fit.model))
                sims[is.nan(sims)] <- 0
                sims[is.infinite(sims)] <- 0
                minsize <- apply(sims,2,cummin)
                QE <- (minsize < Nx)
                ext.boot[,i] <- apply(QE,1,sum)/nsim
            }
        } 


        # Calculate and plot the CI
        ext.ci <- apply(ext.boot,1,quantile,c(alpha/2,1-alpha/2),na.rm=T)
        matplot(t(ext.ci), type='l', lty=2, xlab="Years", ylab="Cumulative extinction risk",col=1,...)
        lines(extrisk)
        output <- list(extrisk=extrisk,year=1:year.max,nboot=nboot,ext.boot=ext.boot,alpha=alpha,ext.ci=t(ext.ci),
                        model.type=paste(model,"count"),model.params=parhat)
        class(output) <- "ext.risk"
        output
    } else if (nboot > 0 && boot.type =="parametric") {
        VC <- summary(fit.model)$cov.unscaled
        if (enoise != "normal") 
            warning("Resampled epsilons and demographic stochasticity not implemented with parametric bootstrap.",
                "Resetting to 'normal'")
        enoise <- "normal"
        boot.pars <- mvrnorm(nboot,coef(fit.model),VC)
        s2 <- var(resid(fit.model))
        ext.boot <- matrix(0,year.max,nboot)
        if (model == "Ricker") {
            boot.pars[,2] <- - boot.pars[,1]/boot.pars[,2]
            for (i in 1:nboot) {
                if (boot.pars[i,2] > 0){
                    sims <- sim.Ricker(Nc, c(boot.pars[i,], s2), year.max, nsim, enoise, resid(model.fit))
                    sims[is.nan(sims)] <- 0
                    sims[is.infinite(sims)] <- 0
                    minsize <- apply(sims,2,cummin)
                    QE <- (minsize < Nx)
                    ext.boot[,i] <- apply(QE,1,sum)/nsim
                }
            }
        } else if (model == "theta") {
            for (i in 1:nboot) {
                if (boot.pars[i,2] > 0){
                    sims <- sim.theta(Nc, c(boot.pars[i,], s2), year.max, nsim, enoise, resid(model.fit))
                    sims[is.nan(sims)] <- 0
                    sims[is.infinite(sims)] <- 0
                    minsize <- apply(sims,2,cummin)
                    QE <- (minsize < Nx)
                    ext.boot[,i] <- apply(QE,1,sum)/nsim
                }
            }
        } else if (model == "logistic") {
            for (i in 1:nboot) {
                if (boot.pars[i,2] > 0){
                    sims <- sim.logistic(Nc, c(boot.pars[i,], s2), year.max, nsim, enoise, resid(model.fit))
                    sims[is.nan(sims)] <- 0
                    sims[is.infinite(sims)] <- 0
                    minsize <- apply(sims,2,cummin)
                    QE <- (minsize < Nx)
                    ext.boot[,i] <- apply(QE,1,sum)/nsim
                }
            }
        } else if (model == "Beverton") {
            for (i in 1:nboot) {
                if (boot.pars[i,2] > 0){
                    sims <- sim.Beverton(Nc, c(boot.pars[i,], s2), year.max, nsim, enoise, resid(model.fit))
                    sims[is.nan(sims)] <- 0
                    sims[is.infinite(sims)] <- 0
                    minsize <- apply(sims,2,cummin)
                    QE <- (minsize < Nx)
                    ext.boot[,i] <- apply(QE,1,sum)/nsim
                }
            }
        } else if (model == "delay") {
            for (i in 1:nboot) {
                if (boot.pars[i,2] < 0){
                    sims <- sim.delay(Nc, c(boot.pars[i,], s2), year.max, nsim, enoise, resid(model.fit))
                    sims[is.nan(sims)] <- 0
                    sims[is.infinite(sims)] <- 0
                    minsize <- apply(sims,2,cummin)
                    QE <- (minsize < Nx)
                    ext.boot[,i] <- apply(QE,1,sum)/nsim
                }
            }
        }

        # Calculate and plot the CI
        ext.ci <- apply(ext.boot,1,quantile,c(alpha/2,1-alpha/2),na.rm=T)
        matplot(t(ext.ci), type='l', lty=2, xlab="Years", ylab="Cumulative extinction risk",ylim=c(0.00000001,max(ext.boot)),col=1,...)
        lines(extrisk)
        output <- list(extrisk=extrisk,year=1:year.max,nboot=nboot,ext.boot=ext.boot,alpha=alpha,ext.ci=t(ext.ci),
                        model.type=paste(model,"count"),model.params=parhat)
        class(output) <- "ext.risk"
        output
    } else {
        output <- list(extrisk=extrisk,year=1:year.max,nboot=NULL,ext.boot=NULL,alpha=NULL,ext.ci=NULL,
                        model.type=paste(model,"count"),model.params=parhat)
        class(output) <- "ext.risk"
        output
    }
}

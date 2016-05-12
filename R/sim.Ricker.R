sim.Ricker <- function(N0, param, Tmax, nsim, method=c("normal","resample"), resids=NULL)
{   
    mu <- param[1]
    K <- param[2]
    if (method[1] == "normal") {
        s <- sqrt(param[3])
        epsilon <- matrix(rnorm(nsim*Tmax,0,s),Tmax,nsim)
    } else if (method[1] == "DS") {
        epsilon <- matrix(rnorm(nsim*Tmax,0,1),Tmax,nsim)
    } else {
        epsilon <- matrix(sample(resids,nsim*Tmax,replace=T),Tmax,nsim)
    }
    Nt <- matrix(N0,1,nsim) 
    for (i in 1:Tmax) {
        if (method[1] == "DS")  {
            s2 <- param[3]+param[4]/Nt[i,]
            s2[s2<0] <- 0
            s2[is.nan(s2)] <- 0
            s2[is.infinite(s2)] <- 0
            s2[Nt[i,] < 0.0000001] <- 0
            epsilon[i,] <- epsilon[i,]*sqrt(s2)
        }
        Nt <- rbind(Nt, Nt[i,]*exp(mu*(1-Nt[i,]/K) + epsilon[i,]))
    }
    Nt[-1,]
}

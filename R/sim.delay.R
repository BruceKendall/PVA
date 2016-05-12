sim.delay <- function(N0, param, Tmax, nsim, method=c("normal","resample"), resids=NULL)
{   
    a0 <- param[1]
    a1 <- param[2]
    a2 <- param[3]
    if (method[1] == "normal") {
        s <- sqrt(param[4])
        epsilon <- matrix(rnorm(nsim*Tmax,0,s),Tmax,nsim)
    } else if (method[1] == "DS") {
        epsilon <- matrix(rnorm(nsim*Tmax,0,1),Tmax,nsim)
    } else {
        epsilon <- matrix(sample(resids,nsim*Tmax,replace=T),Tmax,nsim)
    }
    Nt <- matrix(N0[1],1,nsim)
    Nt <- rbind(Nt,matrix(N0[2],1,nsim))
    for (i in (1:Tmax)+1) {
        if (method[1] == "DS")  {
            s2 <- param[4]+param[5]/Nt[i,]
            s2[s2<0] <- 0
            epsilon[i,] <- epsilon[i,]*sqrt(s2)
        }
        Nt <- rbind(Nt, (Nt[i,]* exp(a0 + a1*Nt[i,] + a2*Nt[i-1,] + epsilon[i-1,])))
    }
    Nt[-c(1,2),]
}

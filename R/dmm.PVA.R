dmm.PVA <- function (x, Nc, Nx, year.max = 100, nsim = 1000, ES = c("MatrixDraw", 
    "ParDraw", "none"), DS = FALSE, Varsamp = c("none", "White"), 
    fecmod = c("logN", "sBeta"), fecmax = 1, ...) 
{
    pbar <- c(x$fec, x$surv, as.vector(x$grow.binom))
    npar <- length(pbar)
    ES <- ES[1]
    Varsamp <- Varsamp[1]
    fecmod = fecmod[1]
    nclass <- length(x$classes)
    if (ES == "MatrixDraw") {
        nmat <- length(x$mat.ann)
        sims <- NULL
        for (sim in 1:nsim) {
            N <- matrix(0, year.max + 1, nclass)
            N[1, ] <- Nc
            imat <- ceiling(runif(year.max, 0, nmat))
            for (yr in 1:year.max) {
                At.dmm <- x$mat.ann[[imat[yr]]]
                if (!DS) {
                  N[yr + 1, ] <- t(At.dmm[[1]] %*% N[yr, ])
                }
                else {
                  if (x$census == "pre") {
                    Pt <- At.dmm$surv[-1]
                    Ft <- At.dmm$fec[-1]
                    P0 <- At.dmm$surv[1]
                  }
                  else {
                    Pt <- At.dmm$surv
                    Ft <- At.dmm$fec
                  }
                  for (i in 1:nclass) {
                    surv <- rbinom(1, N[yr, i], Pt[i])
                    z <- rmultinom(1, surv, At.dmm$grow[, i])
                    if (any(is.na(z))) {
                      s1 <- sum(z, na.rm = T)
                      vec <- is.na(z)
                      z[vec] <- (surv - s1)/sum(vec)
                    }
                    N[yr + 1, ] <- N[yr + 1, ] + z
                    if (x$census == "pre") {
                      fruits <- rpois(1, N[yr, i] * Ft[i])
                      N[yr + 1, 1] <- N[yr + 1, 1] + rbinom(1, 
                        fruits, P0)
                    }
                  }
                  if (x$census == "post") {
                    for (i in 1:nclass) {
                      N[yr + 1, 1] <- N[yr + 1, 1] + rpois(1, 
                        N[yr + 1, i] * Ft[i])
                    }
                  }
                }
            }
            sims <- cbind(sims, apply(N, 1, sum))
        }
    }
    else if (ES == "ParDraw") {
        require(MASS)
        nfec <- length(x$fec)
        nsurv <- length(x$surv)
        nclass <- length(x$classes)
        ngrow <- nclass^2
        if (Varsamp == "none") {
            temp <- VC.dmm(x)
            VC <- temp$VC
            cormat <- temp$cormat
            pbar <- temp$pbar
            Vars <- vector("numeric", npar)
            for (i in 1:npar) Vars[i] <- VC[i, i]
        }
        else {
            temp <- VC.dmm(x)
            cormat <- temp$cormat
            pbar <- temp$pbar
            Vars <- VWhite(x)
            Vars <- c(Vars$VCfec, Vars$VCsurv, Vars$VCgrow)
        }
        Vars[Vars < 1e-12] <- 0
        sims <- NULL
        for (sim in 1:nsim) {
            N <- matrix(0, year.max + 1, nclass)
            N[1, ] <- Nc
            randmat <- mvrnorm(year.max, rep(0, npar), cormat)
            for (i in 1:nfec) {
                if (fecmod == "logN") {
                  mup <- log(pbar[i]) - log(Vars[i]/pbar[i]^2 + 
                    1)/2
                  varp <- log(Vars[i]/pbar[i]^2 + 1)
                  randmat[, i] <- exp(randmat[, i] * sqrt(varp) + 
                    mup)
                }
                else {
                  if (Vars[i] > 0) {
                    p <- pbar[i]/fecmax
                    v <- Vars[i]/fecmax^2
                    temp <- p * (1 - p)/v - 1
                    a <- P * temp
                    b <- (1 - p) * temp
                    randmat[, i] <- qbeta(pnorm(randmat[, i]), 
                      a, b) * fecmax
                  }
                  else {
                    randmat[, i] <- pbar[i]
                  }
                }
            }
            for (i in (nfec + 1):npar) {
                if (Vars[i] > 0) {
                  v <- Vars[i]
                  v <- min(v, pbar[i] * (1 - pbar[i]) - 1e-05)
                  temp <- pbar[i] * (1 - pbar[i])/v - 1
                  a <- pbar[i] * temp
                  b <- (1 - pbar[i]) * temp
                  randmat[, i] <- qbeta(pnorm(randmat[, i]), 
                    a, b)
                }
                else {
                  randmat[, i] <- pbar[i]
                }
            }
            randmat[, pbar == 0] <- 0
            temp <- randmat[, -(1:nfec)]
            temp[, pbar[-(1:nfec)] == 1] <- 1
            randmat[, -(1:nfec)] <- temp
            if (fecmod == "sBeta") {
                temp <- randmat[, 1:nfec]
                temp[, pbar == fecmax] <- fecmax
                randmat[, 1:nfec] <- temp
            }
            for (yr in 1:year.max) {
                if (x$census == "pre") {
                  p0 <- randmat[yr, nfec + 1]
                  lt <- cbind(rep(1, nclass), randmat[yr, (nfec + 
                    2):(nfec + nsurv)], randmat[yr, 2:nfec])
                  gmat <- matrix(randmat[yr, -(1:(nfec + nsurv))], 
                    nclass, nclass)
                }
                else {
                  lt <- cbind(rep(1, nclass), randmat[yr, (nfec + 
                    1):(nsurv + nfec)], randmat[yr, 1:nfec])
                  gmat <- matrix(randmat[yr, -(1:(nfec + nsurv))], 
                    nclass, nclass)
                  p0 <- NULL
                }
                growmat <- matrix(0, nclass, nclass)
                growmat[1, ] <- gmat[1, ]
                for (i in 2:(nclass - 1)) {
                  growmat[i, ] <- gmat[i, ] * (1 - apply(growmat, 
                    2, sum))
                }
                growmat[nclass, ] <- 1 - apply(growmat, 2, sum)
                At.dmm <- dmm(lt, grow = growmat, p0 = p0, classes = x$classes, 
                  census = x$census, data.type = "rates", type = x$type)
                At <- At.dmm$mat
                if (!DS) {
                  N[yr + 1, ] <- t(At %*% N[yr, ])
                }
                else {
                  if (x$census == "pre") {
                    Pt <- At.dmm$surv[-1]
                    Ft <- At.dmm$fec[-1]
                    P0 <- At.dmm$surv[1]
                  }
                  else {
                    Pt <- At.dmm$surv
                    Ft <- At.dmm$fec
                  }
                  for (i in 1:nclass) {
                    surv <- rbinom(1, N[yr, i], Pt[i])
                    z <- rmultinom(1, surv, At.dmm$grow[, i])
                    if (any(is.na(z))) {
                      s1 <- sum(z, na.rm = T)
                      vec <- is.na(z)
                      z[vec] <- (surv - s1)/sum(vec)
                    }
                    N[yr + 1, ] <- N[yr + 1, ] + z
                    if (x$census == "pre") {
                      fruits <- rpois(1, N[yr, i] * Ft[i])
                      N[yr + 1, 1] <- N[yr + 1, 1] + rbinom(1, 
                        fruits, P0)
                    }
                  }
                  if (x$census == "post") {
                    for (i in 1:nclass) {
                      N[yr + 1, 1] <- N[yr + 1, 1] + rpois(1, 
                        N[yr + 1, i] * Ft[i])
                    }
                  }
                }
            }
            sims <- cbind(sims, apply(N, 1, sum))
        }
    }
    else if (ES == "none") {
        if (!DS) 
            stop("Must have either ES or DS")
        At.dmm <- x
        if (x$census == "pre") {
            Pt <- At.dmm$surv[-1]
            Ft <- At.dmm$fec[-1]
            P0 <- At.dmm$surv[1]
        }
        else {
            Pt <- At.dmm$surv
            Ft <- At.dmm$fec
        }
        sims <- NULL
        for (sim in 1:nsim) {
            N <- matrix(0, year.max + 1, nclass)
            N[1, ] <- Nc
            for (yr in 1:year.max) {
                for (i in 1:nclass) {
                  surv <- rbinom(1, N[yr, i], Pt[i])
                  z <- rmultinom(1, surv, At.dmm$grow[, i])
                  if (any(is.na(z))) {
                    s1 <- sum(z, na.rm = T)
                    vec <- is.na(z)
                    z[vec] <- (surv - s1)/sum(vec)
                  }
                  N[yr + 1, ] <- N[yr + 1, ] + z
                  if (x$census == "pre") {
                    fruits <- rpois(1, N[yr, i] * Ft[i])
                    N[yr + 1, 1] <- N[yr + 1, 1] + rbinom(1, 
                      fruits, P0)
                  }
                }
                if (x$census == "post") {
                  for (i in 1:nclass) {
                    N[yr + 1, 1] <- N[yr + 1, 1] + rpois(1, N[yr + 
                      1, i] * Ft[i])
                  }
                }
            }
            sims <- cbind(sims, apply(N, 1, sum))
        }
    }
    sims <- sims[-1, ]
    minsize <- apply(sims, 2, cummin)
    QE <- (minsize < Nx)
    extrisk <- apply(QE, 1, sum)/nsim
    plot(1:year.max, extrisk, type = "l", xlab = "Years", ylab = "Cumulative extinction risk", 
        ...)
    output <- list(extrisk = extrisk, year = 1:year.max, nboot = 0, 
        ext.boot = NULL, alpha = NULL, ext.ci = NULL, model.type = "Demographic", 
        model.params = NULL)
    class(output) <- "ext.risk"
    output
}

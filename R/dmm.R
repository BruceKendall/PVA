dmm <- function (x, fec = NULL, grow = NULL, p0 = NULL, classes = NULL, 
    type = c("age", "stage", "size", "sfa"), census = c("pre", 
        "post"), terminal = FALSE, data.type = c("rates", "aggregated", 
        "raw"), fec.sampling = c("Poisson", "binomial"), seedbank = FALSE, 
    seed.germ = NULL, seed.surv = NULL, annual = F, years = NULL) 
{
    type <- type[1]
    census <- census[1]
    data.type <- data.type[1]
    fec.sampling <- fec.sampling[1]
    if (dim(x)[2] == 3) 
        x <- cbind(rep(1, dim(x)[1]), x)
    if (type == "age") {
        if (data.type == "rates") {
            ages <- sort(unique(x[, 2]))
            nage <- length(ages)
            x1 <- x
            x <- matrix(0, nage, 3)
            for (i in 1:nage) {
                x[i, 1] <- ages[i]
                vec <- x1[, 2] == ages[i]
                x[i, 2] <- mean(x1[vec, 3])
                x[i, 3] <- mean(x1[vec, 4])
            }
            if (census == "post") {
                mat <- matrix(0, nage, nage)
                grow <- matrix(0, nage, nage)
                for (i in 1:(nage - 1)) {
                  mat[i + 1, i] <- x[i, 2]
                  grow[i + 1, i] <- 1
                  mat[1, i] <- x[i, 2] * x[i + 1, 3]
                }
                classes <- as.character(x[, 1])
            }
            else {
                nage <- nage - 1
                mat <- matrix(0, nage, nage)
                grow <- matrix(0, nage, nage)
                for (i in 1:(nage - 1)) {
                  mat[i + 1, i] <- x[i + 1, 2]
                  grow[i + 1, i] <- 1
                  mat[1, i] <- x[1, 2] * x[i + 1, 3]
                }
                mat[1, nage] <- x[1, 2] * x[nage + 1, 3]
                classes <- as.character(x[-1, 1])
            }
            grow[nage, nage] <- 1
            if (!terminal) {
                mat[nage, nage] <- mat[nage, nage - 1]
                if (census == "post") 
                  mat[1, nage] <- x[nage, 2] * x[nage, 3]
                classes[nage] <- paste(classes[nage], "+", sep = "")
            }
            if (annual) {
                mat.ann <- list()
                nyr <- length(years)
                grow <- array(grow, c(nage, nage, nyr))
                for (i in 1:nyr) {
                  vec <- x1[, 1] == years[i]
                  xt <- x1[vec, ]
                  dmmt <- dmm(xt, fec, grow, p0, classes, type, 
                    census, terminal, data.type, fec.sampling, 
                    seedbank, seed.germ, seed.surv)
                  mat.ann <- c(mat.ann, mat = list(dmmt))
                  names(mat.ann)[i] <- paste("dmm", years[i], 
                    sep = "")
                }
            }
            else {
                mat.ann <- NULL
            }
            output <- list(mat = mat, classes = classes, type = type, 
                census = census, surv = x[, 2], fec = x[, 3], 
                grow = grow, grow.binom = NULL, se = NULL, mat.ann = mat.ann)
            class(output) <- "dmm"
            output
        }
    }
    else if (type == "size") {
        if (data.type == "aggregated") {
            from.class <- sort(unique(x[, 2]))
            nclass <- length(from.class)
            pbar <- vector("numeric", nclass)
            gbar <- matrix(0, nclass, nclass)
            fbar <- vector("numeric", nclass)
            N <- vector("numeric", nclass)
            for (i in from.class) {
                ivec <- x[, 2] == i
                jvec <- x[, 3] > 0
                N[i] <- sum(x[ivec, 4])
                Ns <- sum(x[ivec & jvec, 4])
                pbar[i] <- Ns/N[i]
                for (j in from.class) {
                  ijvec <- ivec & x[, 3] == j
                  Ng <- sum(x[ijvec, 4])
                  gbar[j, i] <- Ng/Ns
                }
                ivec <- fec[, 2] == i
                fbar[i] <- sum(fec[ivec, 4])/sum(fec[ivec, 3])
            }
            pbar.se <- sqrt(pbar * (1 - pbar)/N)
            if (fec.sampling == "Poisson") {
                fbar.se <- sqrt(fbar/N)
            }
            else {
                fbar.se <- sqrt(fbar * (1 - fbar)/N)
            }
            gbar.binom <- matrix(0, nclass, nclass)
            gbar.binom.se <- matrix(0, nclass, nclass)
            Nc <- N * pbar
            gbar.binom[1, ] <- gbar[1, ]
            gbar.binom.se[1, ] <- sqrt(gbar.binom[1, ] * (1 - 
                gbar.binom[1, ])/Nc)
            Nc <- Nc * (1 - gbar.binom[1, ])
            vec <- gbar[1, ] < 1
            gbar.binom[2, vec] <- gbar[2, vec]/(1 - gbar[1, vec])
            gbar.binom[2, !vec] <- 1
            gbar.binom[2, gbar.binom[2, ] > 1] <- 1
            j <- Nc > 0
            gbar.binom.se[2, j] <- sqrt(gbar.binom[2, j] * (1 - 
                gbar.binom[2, j])/Nc[j])
            Nc <- Nc * (1 - gbar.binom[2, ])
            for (i in 3:nclass) {
                temp <- apply(gbar[1:(i - 1), ], 2, sum)
                vec <- temp == 1
                temp[vec] <- 0
                gbar.binom[i, ] <- gbar[i, ]/(1 - temp)
                gbar.binom[i, gbar.binom[i, ] > 1] <- 1
                j <- Nc > 0
                gbar.binom.se[i, j] <- sqrt(gbar.binom[i, j] * 
                  (1 - gbar.binom[i, j])/Nc[j])
                Nc <- Nc * (1 - gbar.binom[i, ])
                gbar.binom[i, vec] <- 1
            }
            mat <- matrix(0, nclass, nclass)
            if (seedbank) {
                pbar <- c(seed.surv, pbar)
                pbar.se <- c(0, pbar.se)
                gbar <- rbind(rep(0, nclass), gbar)
                gbar <- cbind(c(1 - seed.germ, seed.germ, rep(0, 
                  nclass - 1)), gbar)
                gbar.binom <- rbind(rep(0, nclass), gbar.binom)
                gbar.binom <- cbind(c(1 - seed.germ, 1, rep(0, 
                  nclass - 1)), gbar.binom)
                gbar.binom.se <- rbind(rep(0, nclass), gbar.binom.se)
                gbar.binom.se <- cbind(rep(0, nclass + 1), gbar.binom.se)
                fbar <- c(0, fbar)
                fbar.se <- c(0, fbar.se)
                nclass <- nclass + 1
                mat <- matrix(0, nclass, nclass)
                for (i in 1:nclass) mat[, i] <- pbar[i] * gbar[, 
                  i]
                if (census == "pre") {
                  mat[1, ] <- mat[1, ] + fbar * seed.surv * (1 - 
                    seed.germ)
                  mat[2, ] <- mat[2, ] + fbar * seed.surv * seed.germ
                }
                else {
                  for (i in 1:nclass) mat[1, i] <- mat[1, i] + 
                    sum(mat[, i] * fbar)
                }
            }
            else {
                for (i in 1:nclass) mat[, i] <- pbar[i] * gbar[, 
                  i]
                if (census == "pre") {
                  mat[1, ] <- mat[1, ] + fbar * p0
                }
                else {
                  for (i in 1:nclass) mat[1, i] <- mat[1, i] + 
                    sum(mat[, i] * fbar)
                }
            }
            if (!is.null(p0)) {
                pbar <- c(p0, pbar)
                fbar <- c(0, fbar)
                pbar.se <- c(0, pbar.se)
                fbar.se <- c(0, fbar.se)
            }
            if (annual) {
                mat.ann <- list()
                nyr <- length(years)
                for (i in 1:nyr) {
                  xt <- x[x[, 1] == years[i], ]
                  ft <- fec[fec[, 1] == years[i], ]
                  dmmt <- dmm(xt, ft, grow, p0, classes, type, 
                    census, terminal, data.type, fec.sampling, 
                    seedbank, seed.germ, seed.surv)
                  mat.ann <- c(mat.ann, mat = list(dmmt))
                  names(mat.ann)[i] <- paste("dmm", years[i], 
                    sep = "")
                }
            }
            else {
                mat.ann <- NULL
            }
            if (is.null(classes)) 
                classes <- as.character(from.class)
            output <- list(mat = mat, classes = classes, type = type, 
                census = census, surv = pbar, fec = fbar, grow = gbar, 
                grow.binom = gbar.binom, se = list(surv = pbar.se, 
                  fec = fbar.se, grow = gbar.binom.se), mat.ann = mat.ann)
            class(output) <- "dmm"
            output
        }
        else if (data.type == "rates") {
            ages <- sort(unique(x[, 2]))
            nage <- length(ages)
            x1 <- x
            x <- matrix(0, nage, 3)
            ng <- dim(grow)[1]
            if (is.matrix(grow)) 
                grow <- array(grow, c(ng, ng, 1))
            g1 <- grow
            grow <- apply(g1, c(1, 2), mean)
            print(g1)
            for (i in 1:nage) {
                x[i, 1] <- ages[i]
                vec <- x1[, 2] == ages[i]
                x[i, 2] <- mean(x1[vec, 3])
                x[i, 3] <- mean(x1[vec, 4])
            }
            pbar <- x[, 2]
            fbar <- x[, 3]
            gbar <- grow
            nclass <- length(pbar)
            if (seedbank) {
                pbar <- c(seed.surv, pbar)
                fbar <- c(0, fbar)
                gbar <- rbind(rep(0, nclass), gbar)
                gbar <- cbind(c(1 - seed.germ, seed.germ, rep(0, 
                  nclass - 1)), gbar)
                nclass <- nclass + 1
            }
            mat <- matrix(0, nclass, nclass)
            for (i in 1:nclass) mat[, i] <- pbar[i] * gbar[, 
                i]
            if (seedbank) {
                if (census == "pre") {
                  mat[1, ] <- mat[1, ] + fbar * seed.surv * (1 - 
                    seed.germ)
                  mat[2, ] <- mat[2, ] + fbar * seed.surv * seed.germ
                }
                else {
                  for (i in 1:nclass) mat[1, i] <- mat[1, i] + 
                    sum(mat[, i] * fbar)
                }
            }
            else {
                for (i in 1:nclass) mat[, i] <- pbar[i] * gbar[, 
                  i]
                if (census == "pre") {
                  mat[1, ] <- mat[1, ] + fbar * p0
                }
                else {
                  for (i in 1:nclass) mat[1, i] <- mat[1, i] + 
                    sum(mat[, i] * fbar)
                }
            }
            if (!is.null(p0)) {
                pbar <- c(p0, pbar)
                fbar <- c(0, fbar)
            }
            if (annual) {
                mat.ann <- list()
                nyr <- length(years)
                for (i in 1:nyr) {
                  xt <- x1[x1[, 1] == years[i], ]
                  growt <- g1[, , i]
                  dmmt <- dmm(xt, fec, growt, p0, classes, type, 
                    census, terminal, data.type, fec.sampling, 
                    seedbank, seed.germ, seed.surv)
                  mat.ann <- c(mat.ann, mat = list(dmmt))
                  names(mat.ann)[i] <- paste("dmm", years[i], 
                    sep = "")
                }
            }
            else {
                mat.ann <- NULL
            }
            if (is.null(classes)) 
                classes <- as.character(from.class)
            output <- list(mat = mat, classes = classes, type = type, 
                census = census, surv = pbar, fec = fbar, grow = gbar, 
                grow.binom = NULL, se = NULL, mat.ann = mat.ann)
            class(output) <- "dmm"
            output
        }
    }
}

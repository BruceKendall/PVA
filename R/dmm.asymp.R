dmm.asymp <- function(x, quiet=F)
{
    if (class(x) == "dmm") {
        mat <- x$mat
        classes <- x$classes
    } else if (class(x) == "matrix") {
        mat <- x
        classes <- 1:(dim(mat)[1])
    } else {
        stop("x must be of class dmm or matrix")
    }
    eigens <- eigen(mat)
    lambda1 <- Re(eigens$value[1])
    stable <- Re(eigens$vector[,1])
    stable <- stable/sum(stable)
    W <- eigens$vector
    V <- Conj(solve(W))
    RV <- Re(V[1,])
    RV <- RV/RV[1]
    sensmat <- Re(V[1,] %*% t(W[,1]))
    sensmat[mat==0] <- 0
    elasmat <- sensmat * mat / lambda1
    
    if (!quiet) {
    cat("Asymptotic matrix analysis\n\n")
    cat("lambda_1 = ",lambda1,"\n\n")
    cat("Class","Stable distribution","Reproductive value\n",sep="\t")
    for (i in 1:length(classes)) 
        cat(classes[i],stable[i],"",RV[i],"\n",sep="\t")
    cat("\nSensitivity matrix:\n")
    print(sensmat)
    cat("\nElasticity matrix:\n")
    print(elasmat)
    }
    invisible(list(lambda1=lambda1, stable=stable, RV=RV, sensmat=sensmat))
}

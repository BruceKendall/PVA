print.dmm <- function(x) {
    nclass <- length(x$classes)
    cat("Demographic ", x$type, "-structured matrix model with ",
        x$census,"-breeding census\n\n",sep="")
    cat("Mean matrix:\n")
    cat("","",x$classes,"\n",sep="\t")
    cat("\t",rep("--------",nclass+1),"\n",sep="")
    for (i in 1:nclass) cat(x$classes[i],"|",round(x$mat[i,],3),"\n",sep="\t")
    cat("\n")
    cat("Mean life table:\n")
    cat("","x","px","fx\n",sep="\t")
    cat("\t",rep("--------",3),"\n",sep="")
    vec <- 1:nclass
    if (x$census=="pre") {
        cat("",0,round(x$surv[1],3),round(x$fec[1],3),"\n",sep="\t")
    }
    for (i in 1:nclass) {
        j <- i
        if (x$census == "pre") j <- j+1
        cat("",x$classes[i],round(x$surv[j],3),round(x$fec[j],3),"\n",sep="\t")
    }
    cat("\n")
    if (x$type != "age") {
        cat("Mean growth matrix:\n")
        cat("","",x$classes,"\n",sep="\t")
        cat("\t",rep("--------",nclass+1),"\n",sep="")
        for (i in 1:nclass) cat(x$classes[i],"|",round(x$grow[i,],3),"\n",sep="\t")
        cat("\n")
    }
}

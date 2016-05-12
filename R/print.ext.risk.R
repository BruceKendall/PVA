print.ext.risk <- function(x, ...) {
    maxrisk <- x$extrisk[length(x$extrisk)]
    cat(x$model.type,"PVA.\n\n")
    cat("Model parameters:\n")
    print(x$model.params)
    cat("\nExtinction risk in", max(x$year), "years is estimated to be", maxrisk,"\n")
    if(!is.null(x$alpha)) cat((1-x$alpha)*100, "% CI is (", x$ext.ci[length(x$extrisk),1],", ", x$ext.ci[length(x$extrisk),2], ")\n",sep="")
    if (maxrisk > 0.5) cat("Median time to extinction is", max(x$year[x$extrisk<0.5]),"years.\n")
    invisible(x)
}

ext.risk <- function(extrisk, year, nboot, ext.boot, alpha, ext.ci, model.type, 
    model.params)
{
    output <- list(extrisk=extrisk, year=year, nboot=nboot, ext.boot=ext.boot, 
        alpha=alpha, ext.ci=ext.ci, model.type=model.type, model.params=model.params)
    class(output) <- "ext.risk"
    output
}

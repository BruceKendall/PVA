se.lambda <- function(x) {
    S <- sens.vr(x)
    sel <- 0
    sel <- sel + sum((S$S.fec * x$se$fec)^2)
    sel <- sel + sum((S$S.surv * x$se$surv)^2)
    if (x$type != "age") 
        sel <- sel + sum((S$S.grow.binom * x$se$grow)^2)
    sqrt(sel)
}

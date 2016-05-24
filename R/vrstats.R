#' Calculate vital rate stats for use in vitalsim
#'
#' @param vrmat A matrix of vital rates, with one year per column; typically the output from \code{\link{extractVRs}}
#'
#' @return A list with elements vrmeans, vrvars, corrin, and corrout, which are used by vitalsim() in the popbio package
#' 
#' @export
#'
vrstats <- function(vrmat) {
  ### Add in some error-checking

  vrmeans <- apply(vrmat, 1, mean)
  vrvars <- apply(vrmat, 1, var)
  corrin <- cor(t(vrmat))
  corrin[is.na(corrin)] <- 0
  corrout <- acf(t(vrmat), lag.max = 1, plot = FALSE)$acf[2, , ]
  corrout[is.nan(corrout)] <- 0
  rownames(corrout) <- rownames(vrmat)
  colnames(corrout) <- rownames(vrmat)

  return(list(vrmeans = vrmeans,
              vrvars = vrvars,
              corrin = corrin,
              corrout = corrout))
}
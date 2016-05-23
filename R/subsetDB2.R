#' Extract subset from a compadre/comadre object
#'
#' @param sub Logical expression that can be used ot subset the metadata component of db
#' @param db A Compadre or Comadre object
#'
#' @details Based on subsetDB() from https://github.com/jonesor/compadreDB; the change is to allow the subsetting to be specified as a logical expression
#' 
#' @return A com(p)adre object (list), with each row of the returned metadata being matched by an element in the other components
#' @export
#'
subsetDB2 <- function(sub, db=comadre){

  #is there a way to put the subset command into the arguments for this function?
#  subsetID <- sub
  e <- substitute(sub)
  r <- eval(e, db$metadata, parent.frame())
  subsetID <- (1:length(r))[r & !is.na(r)]

  # First make a copy of the database.
  ssdb <- db

  # Subset the sub-parts of the database
  ssdb$metadata <- ssdb$metadata[subsetID,]
  ssdb$mat <- ssdb$mat[subsetID]
  ssdb$matrixClass <- ssdb$matrixClass[subsetID]

  # Version information is retained, but modified as follows.
  if("version" %in% names(ssdb)){
    ssdb$version$Version <- paste(ssdb$version$Version," - subset created on ",format(Sys.time(), "%b_%d_%Y"),sep="")
    ssdb$version$DateCreated <- paste(ssdb$version$DateCreated," - subset created on ",format(Sys.time(), "%b_%d_%Y"),sep="")
    ssdb$version$NumberAcceptedSpecies <- length(unique(ssdb$metadata$SpeciesAccepted))
    ssdb$version$NumberStudies <- length(unique(ssdb$metadata$SpeciesAuthor))
    ssdb$version$NumberMatrices <- length(ssdb$mat)
  }

  return(ssdb)
}
#' Reconstitute matrix from vital rates
#'
#' @param vrs A vector or 1-column matrix of vital rates
#' 
#' @details An object called \code{matstruct} (part of the output from \code{\link{extractVRs}}) needs to be in the global environment. This is not passed in because \code{vitalsim} does not allow for the passing of extra parameters through to the matrix construction function.
#'
#' @return A projection matrix
#' @export
vrtomat <- function(vrs) {#}, matstruct) {
  with(matstruct, {
    # This makes available the following variables:
    #   matrank
    #   Findices
    #   Sindices
    #   gtype
    #   Gindices
    #   Gbase

    # There are probably more effient ways of doing this using apply
    # It would be even more efficient if I assumed that the order of the vrs exactly
    #   matches the order in the indices; I'm not sure if that's valid.
    # Even better would be to write a function factory that allows the mapping from vrs to
    #   matrix elements to be created once, rather than at every iteration of the model

    if (!is.vector(vrs)) {
    vrs <- as.vector(vrs)
    names(vrs) <- c(as.character(Findices$element), as.character(Sindices$element), Gindices$element)
    }
 #   print(vrs)
    # Reproduction
    matF <- matrix(0, matrank, matrank)
    Fvec <- grep("F", names(vrs))
    Fvr <- vrs[Fvec]
    for (k in 1:length(Fvr)){
      Findex <- match(names(Fvr)[k], Findices$element)
      matF[Findices$i[Findex], Findices$j[Findex]] <- Fvr[k]
    }

    # Survival and growth
    matU <- matrix(0, matrank, matrank)
    Svec <- grep("S", names(vrs))
    Svr <- vrs[Svec]
    Sj <- Sindices$j
      #as.numeric(stringi::stri_extract_all(Snames, regex="[0-9]", simplify=TRUE))
    for (j in Sj) {
      k <- match(j, Sj)
      if (gtype[j] == "age") { # All survivors progress to next class
        matU[j+1, j] <- Svr[k]
      } else if (gtype[j] == "stasis") { # All survivors stay in same class
        matU[j, j] <- Svr[k]
      } else if (gtype[j] == "stage") { # All stay in same class or advance to next
        jj <- j
        Gindex <- subset(Gindices, j == jj)$element
        g <- vrs[match(names(vrs), Gindex)]
        matU[j, j] <- (1 - g) * Svr[k]
        matU[j+1, j] <- g * Svr[k]
      } else if (gtype[j] == "size") { # Arbitrary growth pattern
        jj <- j
        Glist <- subset(Gindices, j == jj)
        Gvr <- vrs[names(vrs) %in% Glist$element]
        # Make sure things are sorted
        Glist <- Glist[order(Glist$element),]
   #     print(Gvr)
        Gvr <- Gvr[order(names(Gvr))]
        gcol <- rep(0, matrank)
        ng <- length(Gvr)
        # Start at the end
        gcol[Glist$i[ng]] <- Gvr[ng]
        if (ng > 1) {
          for (ii in (ng-1):1) {
            gcol[Glist$i[ii]] <- Gvr[ii] - Gvr[ii+1]
          }
        }
        gcol[Gbase[j]] <- 1 - Gvr[1]
        matU[ , j] <- gcol * Svr[k]
      }
    }
    return(matF + matU)
  })
}